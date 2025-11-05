#!/usr/bin/env python3
"""
Train SAE on Tahoe-100M dataset.

This script trains a large-scale SAE on the Tahoe-100M dataset (95.6M cells)
using the Tahoe X1 3B model at layer 16 with d_sae=40,960 features.
"""

import argparse
import yaml
from pathlib import Path
from typing import Dict, Any
import torch
from torch.utils.data import DataLoader, Dataset, random_split
from rich.console import Console
from rich.progress import Progress
from datasets import load_from_disk
import numpy as np

from sae_genomics.models.tahoe_adapter import TahoeModelAdapter
from sae_genomics.training.trainer import SAETrainer
from sae_genomics.training.config import TrainingConfig

console = Console()


class Tahoe100MDataset(Dataset):
    """Dataset wrapper for Tahoe-100M HuggingFace dataset.

    Converts the HuggingFace format (genes/expressions arrays) to a format
    compatible with Tahoe X1 model.
    """

    def __init__(self, hf_dataset, vocab, max_length: int = 1024):
        """Initialize dataset.

        Args:
            hf_dataset: HuggingFace dataset with 'genes' and 'expressions' columns
            vocab: Gene vocabulary mapping (Ensembl ID -> token ID)
            max_length: Maximum sequence length
        """
        self.dataset = hf_dataset
        self.vocab = vocab
        self.max_length = max_length
        self.vocab_size = len(vocab)  # Store vocab size for validation

    def __len__(self):
        return len(self.dataset)

    def __getitem__(self, idx):
        """Get a single cell.

        Returns:
            Dict with 'gene_ids' and 'expressions' tensors
        """
        sample = self.dataset[idx]

        # Get genes and expressions (these are arrays in the dataset)
        genes = sample['genes']  # Array of gene indices
        expressions = sample['expressions']  # Array of expression values

        # Clip genes to valid vocabulary range (0 to vocab_size-1)
        # Gene IDs >= vocab_size are replaced with -2 (padding token)
        genes = [g if 0 <= g < self.vocab_size else -2 for g in genes[:self.max_length]]
        expressions = expressions[:self.max_length]

        # Convert to torch tensors
        gene_ids = torch.tensor(genes, dtype=torch.long)
        expr_values = torch.tensor(expressions, dtype=torch.float32)

        # Pad if necessary
        if len(gene_ids) < self.max_length:
            pad_len = self.max_length - len(gene_ids)
            gene_ids = torch.cat([gene_ids, torch.full((pad_len,), -2, dtype=torch.long)])
            expr_values = torch.cat([expr_values, torch.zeros(pad_len, dtype=torch.float32)])

        return {
            'gene_ids': gene_ids,
            'expressions': expr_values,
        }


def collate_fn(batch):
    """Collate function for DataLoader."""
    gene_ids = torch.stack([item['gene_ids'] for item in batch])
    expressions = torch.stack([item['expressions'] for item in batch])

    return {
        'gene_ids': gene_ids,
        'expressions': expressions,
    }


def load_config(config_path: Path) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description="Train SAE on Tahoe-100M dataset")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/training/tahoe100m_3b_layer16.yaml"),
        help="Path to training configuration file",
    )
    parser.add_argument(
        "--data-path",
        type=Path,
        default=None,
        help="Override dataset path from config",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Override output directory from config",
    )
    parser.add_argument(
        "--resume",
        type=Path,
        default=None,
        help="Resume from checkpoint",
    )

    args = parser.parse_args()

    # Load configuration
    console.print("[bold green]Loading configuration...[/bold green]")
    config = load_config(args.config)

    # Override config with CLI arguments
    if args.data_path:
        config['data_path'] = str(args.data_path)
    if args.output_dir:
        config['checkpoint_dir'] = str(args.output_dir)

    console.print(f"Configuration: {args.config}")
    console.print(f"Model: Tahoe X1 {config['model_size']} (d_model={2560 if config['model_size']=='3b' else 512})")
    console.print(f"Layer: {config['layer']}")
    console.print(f"Hook point: {config['hook_point']}")
    console.print(f"SAE dimensions: d_in={2560 if config['model_size']=='3b' else 512}, d_sae={config['d_sae']}, k={config['k']}")
    console.print(f"Training steps: {config['steps']}")
    console.print(f"Output directory: {config['checkpoint_dir']}\n")

    # Create output directory
    output_dir = Path(config['checkpoint_dir'])
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Load Tahoe-100M dataset
    console.print("[bold]Step 1: Loading Tahoe-100M dataset[/bold]")
    data_path = Path(config['data_path'])

    if not data_path.exists():
        console.print(f"[red]Error: Dataset not found at {data_path}[/red]")
        console.print(f"[yellow]Please download the dataset first using:[/yellow]")
        console.print(f"  python scripts/download_tahoe100m.py")
        return 1

    console.print(f"Loading dataset from {data_path}...")
    hf_dataset = load_from_disk(str(data_path))
    console.print(f"✓ Loaded dataset with {len(hf_dataset):,} cells")
    console.print(f"  Columns: {hf_dataset.column_names}")
    console.print(f"  Features: {hf_dataset.features}\n")

    # Step 2: Load Tahoe X1 model
    console.print("[bold]Step 2: Loading Tahoe X1 model[/bold]")
    tahoe_adapter = TahoeModelAdapter.from_hf(
        model_size=config['model_size'],
        device=config['device'],
    )
    console.print(f"✓ Loaded {config['model_size']} model")
    console.print(f"  Vocabulary size: {len(tahoe_adapter.vocab):,}\n")

    # Step 3: Create train/val split
    console.print("[bold]Step 3: Creating train/validation split[/bold]")
    train_size = int(len(hf_dataset) * config['train_split'])
    val_size = len(hf_dataset) - train_size

    console.print(f"Splitting dataset: {config['train_split']:.0%} train, {config['val_split']:.0%} val")
    train_dataset_hf, val_dataset_hf = random_split(
        hf_dataset,
        [train_size, val_size],
        generator=torch.Generator().manual_seed(42),
    )

    console.print(f"✓ Train: {len(train_dataset_hf):,} cells")
    console.print(f"✓ Val: {len(val_dataset_hf):,} cells\n")

    # Wrap in custom dataset
    train_dataset = Tahoe100MDataset(
        train_dataset_hf,
        tahoe_adapter.vocab,
        max_length=config['max_length'],
    )
    val_dataset = Tahoe100MDataset(
        val_dataset_hf,
        tahoe_adapter.vocab,
        max_length=config['max_length'],
    )

    # Create dataloaders
    train_loader = DataLoader(
        train_dataset,
        batch_size=config['batch_size'],
        shuffle=True,
        num_workers=config['num_workers'],
        collate_fn=collate_fn,
        pin_memory=True,
    )
    val_loader = DataLoader(
        val_dataset,
        batch_size=config['batch_size'],
        shuffle=False,
        num_workers=config['num_workers'],
        collate_fn=collate_fn,
        pin_memory=True,
    )

    console.print(f"✓ Created dataloaders")
    console.print(f"  Train batches: {len(train_loader):,}")
    console.print(f"  Val batches: {len(val_loader):,}\n")

    # Step 4: Initialize SAE trainer
    console.print("[bold]Step 4: Initializing SAE trainer[/bold]")

    # Determine d_in from model size
    d_in = 2560 if config['model_size'] == '3b' else (2048 if config['model_size'] == '1.3b' else 512)

    # Create training config
    training_config = TrainingConfig(
        lr=config['lr'],
        batch_size=config['batch_size'],
        total_training_steps=config['steps'],
        eval_every=config['eval_interval'],
        save_every=config['save_interval'],
        log_every=config['log_interval'],
        device=config['device'],
        d_in=d_in,
        d_sae=config['d_sae'],
        k=config['k'],
        activation=config['activation'],
    )

    # Initialize trainer
    if args.resume:
        console.print(f"Resuming from checkpoint: {args.resume}")
        trainer = SAETrainer.load_checkpoint(args.resume, device=config['device'])
    else:
        trainer = SAETrainer(
            d_in=d_in,
            d_sae=config['d_sae'],
            k=config['k'],
            activation=config['activation'],
            lr=config['lr'],
            l1_coefficient=0.0,  # TopK doesn't need L1
            device=config['device'],
        )

    console.print(f"✓ Initialized SAE trainer")
    console.print(f"  d_in: {trainer.d_in}")
    console.print(f"  d_sae: {trainer.d_sae}")
    console.print(f"  k: {trainer.k}")
    console.print(f"  activation: {trainer.activation}\n")

    # Step 5: Extract activations from Tahoe X1
    console.print("[bold]Step 5: Extracting activations from Tahoe X1[/bold]")
    console.print(f"Extracting from layer {config['layer']} at {config['hook_point']}")
    console.print("This will take a while for 95.6M cells...\n")

    # Set up hook point for layer 16 MLP output
    hook_point = f"model.transformer_encoder.layers.{config['layer']}"
    tahoe_adapter.hook_points = [hook_point]

    # Extract activations in batches
    console.print("Processing training set...")
    train_activations_list = []
    train_gene_ids_list = []

    with Progress() as progress:
        task = progress.add_task("[cyan]Extracting train activations...", total=len(train_loader))

        for batch_idx, batch in enumerate(train_loader):
            # Register hooks
            tahoe_adapter.register_hooks()

            # Move batch to device
            batch = {k: v.to(tahoe_adapter.device) if isinstance(v, torch.Tensor) else v
                    for k, v in batch.items()}

            # Convert to bfloat16 for flash attention
            if tahoe_adapter.device.type == 'cuda':
                batch = {k: v.to(torch.bfloat16) if isinstance(v, torch.Tensor) and v.dtype == torch.float32 else v
                        for k, v in batch.items()}

            # Rename keys to match Tahoe expectations
            # Create gen_mask: True for valid genes, False for padding
            gen_mask = batch['gene_ids'] != -2  # -2 is padding token
            tahoe_batch = {
                'gene': batch['gene_ids'],
                'expr': batch['expressions'],
                'gen_mask': gen_mask,
            }

            # Forward pass to extract activations
            with torch.no_grad():
                _ = tahoe_adapter.model(tahoe_batch)

            # Get activations
            if hook_point in tahoe_adapter.activations:
                act = tahoe_adapter.activations[hook_point].cpu()
                if act.dtype == torch.bfloat16:
                    act = act.to(torch.float32)
                train_activations_list.append(act)
                train_gene_ids_list.append(batch['gene_ids'].cpu())

            # Clear for next batch
            tahoe_adapter.activations = {}
            tahoe_adapter.remove_hooks()

            progress.update(task, advance=1)

            # Limit for memory constraints (full 95.6M cells would require ~petabytes)
            # Using 10,000 batches = 320K cells = ~320GB of activations
            if batch_idx >= 10000:  # 320,000 cells (0.3% of full dataset)
                console.print(f"[cyan]Processed {batch_idx+1} batches ({(batch_idx+1)*32:,} cells) for training[/cyan]")
                break

    # Concatenate all batches
    console.print("Concatenating activations...")
    train_activations = torch.cat(train_activations_list, dim=0)
    train_gene_ids = torch.cat(train_gene_ids_list, dim=0)

    console.print(f"✓ Extracted train activations: {train_activations.shape}")
    console.print(f"  Shape: (cells={train_activations.shape[0]}, seq_len={train_activations.shape[1]}, d_model={train_activations.shape[2]})\n")

    # Extract validation set
    console.print("Processing validation set...")
    val_activations_list = []
    val_gene_ids_list = []

    with Progress() as progress:
        task = progress.add_task("[cyan]Extracting val activations...", total=min(len(val_loader), 1000))

        for batch_idx, batch in enumerate(val_loader):
            if batch_idx >= 1000:  # 32,000 validation cells
                break

            tahoe_adapter.register_hooks()
            batch = {k: v.to(tahoe_adapter.device) if isinstance(v, torch.Tensor) else v
                    for k, v in batch.items()}
            if tahoe_adapter.device.type == 'cuda':
                batch = {k: v.to(torch.bfloat16) if isinstance(v, torch.Tensor) and v.dtype == torch.float32 else v
                        for k, v in batch.items()}

            gen_mask = batch['gene_ids'] != -2  # -2 is padding token
            tahoe_batch = {'gene': batch['gene_ids'], 'expr': batch['expressions'], 'gen_mask': gen_mask}

            with torch.no_grad():
                _ = tahoe_adapter.model(tahoe_batch)

            if hook_point in tahoe_adapter.activations:
                act = tahoe_adapter.activations[hook_point].cpu()
                if act.dtype == torch.bfloat16:
                    act = act.to(torch.float32)
                val_activations_list.append(act)
                val_gene_ids_list.append(batch['gene_ids'].cpu())

            tahoe_adapter.activations = {}
            tahoe_adapter.remove_hooks()
            progress.update(task, advance=1)

    val_activations = torch.cat(val_activations_list, dim=0)
    val_gene_ids = torch.cat(val_gene_ids_list, dim=0)

    console.print(f"✓ Extracted val activations: {val_activations.shape}\n")

    # Save activations to disk
    console.print("Saving activations to disk...")
    torch.save({
        'train_activations': train_activations,
        'train_gene_ids': train_gene_ids,
        'val_activations': val_activations,
        'val_gene_ids': val_gene_ids,
    }, output_dir / 'activations.pt')
    console.print(f"✓ Saved to {output_dir / 'activations.pt'}\n")

    # Step 6: Train SAE
    console.print("[bold]Step 6: Training SAE[/bold]")
    console.print(f"Training for {config['steps']:,} steps...")
    console.print(f"Checkpoints will be saved to: {output_dir}\n")

    # Flatten activations for SAE training
    # Shape: (cells, seq_len, d_model) -> (cells * seq_len, d_model)
    train_acts_flat = train_activations.reshape(-1, train_activations.shape[-1])
    val_acts_flat = val_activations.reshape(-1, val_activations.shape[-1])

    console.print(f"Flattened activations:")
    console.print(f"  Train: {train_acts_flat.shape}")
    console.print(f"  Val: {val_acts_flat.shape}\n")

    # Train SAE using the trainer
    trainer.train_on_activations(
        activations=train_acts_flat,
        n_epochs=config['steps'] // (len(train_acts_flat) // config['batch_size']),
        batch_size=config['batch_size'],
        log_every=config['log_interval'],
        save_every=config['save_interval'],
        output_dir=output_dir,
    )

    console.print("\n[bold green]✓ Training complete![/bold green]")
    console.print(f"Results saved to: {output_dir}")
    console.print(f"  - SAE checkpoint: {output_dir / 'sae_final.pt'}")
    console.print(f"  - Training metrics: {output_dir / 'training_metrics.json'}")
    console.print(f"  - Activations: {output_dir / 'activations.pt'}")

    return 0


if __name__ == "__main__":
    exit(main())
