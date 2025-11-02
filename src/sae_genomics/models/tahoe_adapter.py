"""Adapter for loading and extracting activations from Tahoe X1 models."""

import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader

# Add tahoe-x1 to path
TAHOE_PATH = Path(__file__).parent.parent.parent.parent / "external" / "tahoe-x1"
if str(TAHOE_PATH) not in sys.path:
    sys.path.insert(0, str(TAHOE_PATH))

from tahoe_x1.model import ComposerTX
from tahoe_x1.utils.util import loader_from_adata


class TahoeModelAdapter:
    """Adapter for Tahoe X1 models to extract intermediate activations."""

    def __init__(
        self,
        model: ComposerTX,
        vocab: Dict[str, int],
        hook_points: Optional[List[str]] = None,
        device: str = "auto",
    ):
        """Initialize the adapter.

        Args:
            model: Loaded Tahoe X1 ComposerTX model
            vocab: Gene vocabulary mapping
            hook_points: List of layer names to extract activations from
            device: Device to run on ('auto', 'cuda', 'cpu')
        """
        self.model = model
        self.vocab = vocab
        self.hook_points = hook_points or self._get_default_hook_points()

        if device == "auto":
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        else:
            self.device = torch.device(device)

        self.model = self.model.to(self.device)
        self.model.eval()

        # Storage for activations
        self.activations = {}
        self.hooks = []

    def _get_default_hook_points(self) -> List[str]:
        """Get default hook points based on model size."""
        n_layers = self.model.model.n_layers

        if n_layers <= 12:  # 70M model
            return [f"model.blocks.{i}" for i in [0, 3, 6, 9, 11]]
        elif n_layers <= 24:  # 1.3B model
            return [f"model.blocks.{i}" for i in [0, 6, 12, 18, 23]]
        else:  # 3B model
            return [f"model.blocks.{i}" for i in [0, 8, 16, 24, 31]]

    def _get_activation_hook(self, name: str):
        """Create a hook function to capture activations."""
        def hook(module, input, output):
            # Store the output activation
            if isinstance(output, tuple):
                self.activations[name] = output[0].detach()
            else:
                self.activations[name] = output.detach()
        return hook

    def register_hooks(self):
        """Register forward hooks to capture activations."""
        self.remove_hooks()  # Clear any existing hooks

        for hook_point in self.hook_points:
            # Navigate to the module
            module = self.model
            for part in hook_point.split('.'):
                module = getattr(module, part)

            # Register hook
            hook = module.register_forward_hook(self._get_activation_hook(hook_point))
            self.hooks.append(hook)

    def remove_hooks(self):
        """Remove all registered hooks."""
        for hook in self.hooks:
            hook.remove()
        self.hooks = []
        self.activations = {}

    @torch.no_grad()
    def extract_activations(
        self,
        dataloader: DataLoader,
        layer: Optional[str] = None,
        max_batches: Optional[int] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Extract activations from the model.

        Args:
            dataloader: DataLoader with single-cell data
            layer: Specific layer to extract from (default: last hook point)
            max_batches: Maximum number of batches to process

        Returns:
            activations: Tensor of shape (n_samples, seq_len, d_model)
            gene_ids: Tensor of shape (n_samples, seq_len) with gene IDs
        """
        if layer is None:
            layer = self.hook_points[-1]  # Use last layer by default

        self.register_hooks()

        all_activations = []
        all_gene_ids = []

        for batch_idx, batch in enumerate(dataloader):
            if max_batches is not None and batch_idx >= max_batches:
                break

            # Move batch to device
            if isinstance(batch, dict):
                batch = {k: v.to(self.device) if isinstance(v, torch.Tensor) else v
                        for k, v in batch.items()}
                input_ids = batch.get('input_ids', batch.get('gene_ids'))
            else:
                batch = batch.to(self.device)
                input_ids = batch

            # Forward pass
            _ = self.model(batch)

            # Extract activation for specified layer
            if layer in self.activations:
                act = self.activations[layer].cpu()
                all_activations.append(act)
                all_gene_ids.append(input_ids.cpu())

            # Clear activations for next batch
            self.activations = {}

        self.remove_hooks()

        # Concatenate all batches
        activations = torch.cat(all_activations, dim=0)
        gene_ids = torch.cat(all_gene_ids, dim=0)

        return activations, gene_ids

    @classmethod
    def from_hf(
        cls,
        model_size: str = "70m",
        repo_id: str = "tahoebio/Tahoe-x1",
        hook_points: Optional[List[str]] = None,
        device: str = "auto",
    ) -> "TahoeModelAdapter":
        """Load model from HuggingFace.

        Args:
            model_size: Model size ('70m', '1.3b', '3b')
            repo_id: HuggingFace repository ID
            hook_points: Specific hook points to use
            device: Device to run on

        Returns:
            TahoeModelAdapter instance
        """
        model, vocab, _, _ = ComposerTX.from_hf(repo_id, model_size)
        return cls(model, vocab, hook_points, device)

    def create_dataloader_from_adata(
        self,
        adata,
        batch_size: int = 32,
        max_length: int = 2048,
        gene_id_key: str = "ensembl_id",
        num_workers: int = 4,
        prefetch_factor: int = 2,
    ) -> DataLoader:
        """Create a DataLoader from AnnData object.

        Args:
            adata: AnnData object with single-cell data
            batch_size: Batch size
            max_length: Maximum sequence length
            gene_id_key: Key in adata.var for gene IDs
            num_workers: Number of data loading workers
            prefetch_factor: Data prefetching factor

        Returns:
            DataLoader
        """
        # Match genes to vocabulary
        adata.var["id_in_vocab"] = [
            self.vocab[gene] if gene in self.vocab else -1
            for gene in adata.var[gene_id_key]
        ]

        gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
        print(f"Matched {np.sum(gene_ids_in_vocab >= 0)}/{len(gene_ids_in_vocab)} genes to vocabulary")

        # Filter to genes in vocabulary
        adata = adata[:, adata.var["id_in_vocab"] >= 0]
        genes = adata.var[gene_id_key].tolist()
        gene_ids = np.array([self.vocab[gene] for gene in genes], dtype=int)

        # Get collator config from model
        collator_cfg = {
            "pad_token_id": getattr(self.model.model, "pad_token_id", 0),
            "mask_token_id": getattr(self.model.model, "mask_token_id", 1),
        }

        # Create DataLoader
        loader = loader_from_adata(
            adata=adata,
            collator_cfg=collator_cfg,
            vocab=self.vocab,
            batch_size=batch_size,
            max_length=max_length,
            gene_ids=gene_ids,
            num_workers=num_workers,
            prefetch_factor=prefetch_factor,
        )

        return loader

    def get_model_config(self) -> Dict:
        """Get model configuration."""
        return {
            "n_layers": self.model.model.n_layers,
            "d_model": self.model.model.d_model,
            "n_heads": self.model.model.n_heads,
            "vocab_size": len(self.vocab),
            "device": str(self.device),
            "hook_points": self.hook_points,
        }
