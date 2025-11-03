"""Adapter for loading and extracting activations from Tahoe X1 models."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader

if TYPE_CHECKING:
    from tahoe_x1.model import ComposerTX

# Add tahoe-x1 to path
TAHOE_PATH = Path(__file__).parent.parent.parent.parent / "external" / "tahoe-x1"
if str(TAHOE_PATH) not in sys.path:
    sys.path.insert(0, str(TAHOE_PATH))

# Lazy imports - only load when needed
def _import_tahoe():
    """Lazy import of tahoe_x1 to avoid dependency issues."""
    from tahoe_x1.model import ComposerTX
    from tahoe_x1.utils.util import loader_from_adata
    return ComposerTX, loader_from_adata


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

        # Convert to bfloat16 for flash attention compatibility (flash attention requires fp16/bf16)
        if self.device.type == 'cuda':
            self.model = self.model.to(torch.bfloat16)

        self.model.eval()

        # Storage for activations
        self.activations = {}
        self.hooks = []

    def _get_default_hook_points(self) -> List[str]:
        """Get default hook points based on model size."""
        n_layers = self.model.model.n_layers

        # Hook into layer outputs (entire block)
        # Tahoe uses: model.transformer_encoder.layers[i]
        if n_layers <= 12:  # 70M model
            return [f"model.transformer_encoder.layers.{i}" for i in [0, 3, 6, 9, 11]]
        elif n_layers <= 24:  # 1.3B model
            return [f"model.transformer_encoder.layers.{i}" for i in [0, 6, 12, 18, 23]]
        else:  # 3B model
            return [f"model.transformer_encoder.layers.{i}" for i in [0, 8, 16, 24, 31]]

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
            # Navigate to the module with proper error handling
            module = self.model
            try:
                for part in hook_point.split('.'):
                    module = getattr(module, part)
            except AttributeError as e:
                # Provide helpful error message with correct structure
                n_layers = self.model.model.n_layers
                available_layers = [f"model.transformer_encoder.layers.{i}" for i in range(n_layers)]
                raise ValueError(
                    f"Invalid hook point '{hook_point}': {e}. "
                    f"Check that this layer exists in the model. "
                    f"Tahoe model structure: model.transformer_encoder.layers[0-{n_layers-1}]. "
                    f"Example hook points: {available_layers[:3]}"
                ) from e

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

            # Move batch to device and extract gene IDs
            if isinstance(batch, dict):
                batch = {k: v.to(self.device) if isinstance(v, torch.Tensor) else v
                        for k, v in batch.items()}
                # Convert to bfloat16 for flash attention compatibility
                if self.device.type == 'cuda':
                    batch = {k: v.to(torch.bfloat16) if isinstance(v, torch.Tensor) and v.dtype == torch.float32 else v
                            for k, v in batch.items()}
                # Tahoe dataloader returns: gene, expr, expr_target, expr_raw, gen_mask
                gene_ids_batch = batch.get('gene', batch.get('input_ids', batch.get('gene_ids')))
                if gene_ids_batch is None:
                    raise ValueError(
                        "Batch dict must contain 'gene', 'input_ids', or 'gene_ids' key. "
                        f"Available keys: {list(batch.keys())}"
                    )
            else:
                if not isinstance(batch, torch.Tensor):
                    raise TypeError(
                        f"Expected batch to be dict or Tensor, got {type(batch).__name__}"
                    )
                batch = {k: v.to(self.device) if isinstance(v, torch.Tensor) else v
                        for k, v in batch.items()}
                if self.device.type == 'cuda':
                    batch = {k: v.to(torch.bfloat16) if isinstance(v, torch.Tensor) and v.dtype == torch.float32 else v
                            for k, v in batch.items()}
                gene_ids_batch = batch

            # Forward pass - Tahoe model expects a dict with gene, expr, etc.
            _ = self.model(batch)

            # Extract activation for specified layer
            if layer in self.activations:
                act = self.activations[layer].cpu()
                # Convert back to float32 for compatibility with SAE training
                if act.dtype == torch.bfloat16:
                    act = act.to(torch.float32)
                all_activations.append(act)
                all_gene_ids.append(gene_ids_batch.cpu())

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
        ComposerTX, _ = _import_tahoe()
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
        _, loader_from_adata = _import_tahoe()

        # Match genes to vocabulary
        adata.var["id_in_vocab"] = [
            self.vocab[gene] if gene in self.vocab else -1
            for gene in adata.var[gene_id_key]
        ]

        gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
        n_matched = np.sum(gene_ids_in_vocab >= 0)
        import sys
        print(f"Matched {n_matched}/{len(gene_ids_in_vocab)} genes to vocabulary", file=sys.stderr)

        # Validate that at least some genes matched
        if n_matched == 0:
            raise ValueError(
                f"No genes in the dataset matched the model vocabulary. "
                f"Check that gene_id_key='{gene_id_key}' contains the correct gene identifiers. "
                f"Example genes in data: {list(adata.var[gene_id_key].head())}"
            )

        # Filter to genes in vocabulary
        adata = adata[:, adata.var["id_in_vocab"] >= 0]
        genes = adata.var[gene_id_key].tolist()
        gene_ids = np.array([self.vocab[gene] for gene in genes], dtype=int)

        # Get collator config from model
        # These values are needed by tahoe's loader_from_adata
        # Must be an OmegaConf DictConfig (not regular dict) for attribute access
        from omegaconf import DictConfig

        collator_cfg = DictConfig({
            # Required fields (accessed with attribute access)
            "pad_token_id": getattr(self.model.model, "pad_token_id", 0),
            "pad_value": getattr(self.model.model, "pad_value", -2),  # Expression value for PAD tokens
            "mask_value": -1,  # Expression value for MASK tokens (not used in inference)
            "mlm_probability": 0.0,  # No masking for inference
            "sampling": True,  # Sample genes when max_length < n_genes
            # Optional fields with sensible defaults for inference
            "do_padding": True,
            "do_binning": True,
            "log_transform": False,
            "num_bins": 51,
            "right_binning": False,
            "keep_first_n_tokens": 1,
            "use_chem_token": False,
        })

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
