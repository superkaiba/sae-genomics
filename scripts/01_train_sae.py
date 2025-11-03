#!/usr/bin/env python3
"""
Train a Sparse Autoencoder on Tahoe X1 activations.

Usage:
    python scripts/01_train_sae.py \\
        --data path/to/data.h5ad \\
        --model-size 70m \\
        --output results/my_experiment \\
        --steps 10000

Example with config:
    python scripts/01_train_sae.py \\
        --data data/example.h5ad \\
        --config configs/training/default.yaml \\
        --output results/exp1
"""

import argparse
import json
from pathlib import Path

import scanpy as sc
import torch
from rich.console import Console
from rich.progress import track

from sae_genomics.models.model_config import ModelConfig
from sae_genomics.models.tahoe_adapter import TahoeModelAdapter
from sae_genomics.training.config import TrainingConfig
from sae_genomics.training.trainer import SAETrainer

console = Console()


def main():
    parser = argparse.ArgumentParser(description="Train SAE on Tahoe X1 activations")
    parser.add_argument(
        "--data",
        type=Path,
        required=True,
        help="Path to AnnData file (.h5ad) with single-cell data",
    )
    parser.add_argument(
        "--model-size",
        type=str,
        default="70m",
        choices=["70m", "1.3b", "3b"],
        help="Tahoe X1 model size",
    )
    parser.add_argument(
        "--model-config",
        type=Path,
        default=None,
        help="Path to model config (overrides --model-size)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Path to training configuration YAML",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default="./results/sae_training",
        help="Output directory for results",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=None,
        help="Training steps (overrides config)",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=None,
        help="Maximum number of cells to use",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="auto",
        choices=["auto", "cuda", "cpu"],
        help="Device to use",
    )
    parser.add_argument(
        "--gene-id-key",
        type=str,
        default="ensembl_id",
        help="Key in adata.var for gene IDs",
    )

    args = parser.parse_args()

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    console.print("[bold green]SAE Training Pipeline[/bold green]")
    console.print(f"Data: {args.data}")
    console.print(f"Model: Tahoe X1 {args.model_size}")
    console.print(f"Output: {args.output}")

    # Load configurations
    if args.config:
        console.print(f"Loading training config from {args.config}")
        train_config = TrainingConfig.from_yaml(args.config)
    else:
        console.print("Using default training configuration")
        train_config = TrainingConfig()

    if args.model_config:
        console.print(f"Loading model config from {args.model_config}")
        model_config = ModelConfig.from_yaml(args.model_config)
    else:
        # Use default config for model size
        default_configs = {
            "70m": "configs/models/tx1_70m.yaml",
            "1.3b": "configs/models/tx1_3b.yaml",  # Both use same config for now
            "3b": "configs/models/tx1_3b.yaml",    # TODO: Create separate configs if needed
        }
        config_path = Path(default_configs[args.model_size])
        if config_path.exists():
            model_config = ModelConfig.from_yaml(config_path)
        else:
            console.print(f"[yellow]Warning: Config file {config_path} not found, using defaults[/yellow]")
            model_config = None

    # Override with command-line arguments
    if args.steps:
        train_config.total_training_steps = args.steps
    if args.max_cells:
        train_config.max_cells = args.max_cells
    if args.device != "auto":
        train_config.device = args.device

    # Save configurations
    with open(args.output / "train_config.json", "w") as f:
        json.dump(train_config.to_dict(), f, indent=2)
    if model_config:
        with open(args.output / "model_config.json", "w") as f:
            json.dump(model_config.to_dict(), f, indent=2)

    # Step 1: Load Tahoe X1 model
    console.print("\n[bold]Step 1: Loading Tahoe X1 model[/bold]")
    tahoe_adapter = TahoeModelAdapter.from_hf(
        model_size=args.model_size,
        device=train_config.device,
    )
    console.print(f"✓ Loaded {args.model_size} model")
    console.print(f"  - Layers: {tahoe_adapter.model.model.n_layers}")
    console.print(f"  - d_model: {tahoe_adapter.model.model.d_model}")
    console.print(f"  - Vocab size: {len(tahoe_adapter.vocab)}")

    # Step 2: Load data
    console.print("\n[bold]Step 2: Loading single-cell data[/bold]")
    adata = sc.read_h5ad(args.data)
    console.print(f"✓ Loaded {adata.n_obs} cells × {adata.n_vars} genes")

    # Subset if requested
    if train_config.max_cells and adata.n_obs > train_config.max_cells:
        console.print(f"  Subsampling to {train_config.max_cells} cells")
        adata = adata[:train_config.max_cells, :]

    # Step 3: Create DataLoader
    console.print("\n[bold]Step 3: Creating DataLoader[/bold]")
    dataloader = tahoe_adapter.create_dataloader_from_adata(
        adata,
        batch_size=train_config.batch_size,
        max_length=model_config.context_length if model_config else 1024,
        gene_id_key=args.gene_id_key,
        num_workers=4,
        prefetch_factor=2,
    )
    console.print(f"✓ Created DataLoader with batch size {train_config.batch_size}")

    # Step 4: Extract activations
    console.print("\n[bold]Step 4: Extracting activations from Tahoe X1[/bold]")
    console.print("This may take a while...")

    # Extract a batch to get d_model
    sample_acts, sample_genes = tahoe_adapter.extract_activations(
        dataloader,
        max_batches=1,
    )
    d_model = sample_acts.shape[-1]
    console.print(f"✓ Activation dimension: {d_model}")

    # Now extract full dataset
    console.print("Extracting full activations...")
    max_batches = min(
        len(dataloader),
        (train_config.total_training_steps * train_config.batch_size) // train_config.batch_size,
    )
    activations, gene_ids = tahoe_adapter.extract_activations(
        dataloader,
        max_batches=max_batches,
    )
    console.print(f"✓ Extracted activations: {activations.shape}")

    # Save activations for later use
    activations_path = args.output / "activations.pt"
    torch.save(
        {
            "activations": activations,
            "gene_ids": gene_ids,
        },
        activations_path,
    )
    console.print(f"✓ Saved activations to {activations_path}")

    # Step 5: Initialize and train SAE
    console.print("\n[bold]Step 5: Training Sparse Autoencoder[/bold]")

    # Update d_in if we got it from model
    if train_config.d_in != d_model:
        console.print(f"  Updating d_in: {train_config.d_in} → {d_model}")
        train_config.d_in = d_model

    trainer = SAETrainer(
        d_in=train_config.d_in,
        d_sae=train_config.d_sae,
        activation=train_config.activation,
        k=train_config.k,
        lr=train_config.lr,
        l1_coefficient=train_config.l1_coefficient,
        device=train_config.device,
    )

    console.print(f"✓ Initialized SAE")
    console.print(f"  - d_in: {train_config.d_in}")
    console.print(f"  - d_sae: {train_config.d_sae}")
    console.print(f"  - Activation: {train_config.activation} (k={train_config.k})")

    # Calculate number of epochs
    sae_batch_size = 256  # SAE training batch size (independent of dataloader batch size)
    n_samples = activations.shape[0] * activations.shape[1]  # cells * seq_len
    steps_per_epoch = n_samples // sae_batch_size
    n_epochs = max(1, train_config.total_training_steps // steps_per_epoch)

    console.print(f"  - Training for {n_epochs} epochs ({train_config.total_training_steps} steps)")
    console.print(f"  - Steps per epoch: {steps_per_epoch} (batch_size={sae_batch_size})")

    # Train
    metrics = trainer.train_on_activations(
        activations=activations,
        n_epochs=n_epochs,
        batch_size=sae_batch_size,  # SAE training batch size
        log_every=train_config.log_every,
        save_every=train_config.save_every,
        output_dir=args.output / "checkpoints",
    )

    console.print(f"✓ Training complete!")

    # Save final checkpoint
    final_checkpoint = args.output / "sae_final.pt"
    trainer.save_checkpoint(final_checkpoint)
    console.print(f"✓ Saved final model to {final_checkpoint}")

    # Save training metrics
    with open(args.output / "training_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)

    # Summary
    console.print("\n[bold green]Training Complete![/bold green]")
    console.print(f"Results saved to: {args.output}")
    console.print(f"Final checkpoint: {final_checkpoint}")
    console.print(f"\nNext steps:")
    console.print(f"  1. Extract feature activations:")
    console.print(f"     python scripts/02_extract_features.py --sae {final_checkpoint} --data {args.data}")
    console.print(f"  2. Create dashboard:")
    console.print(f"     python scripts/03_create_dashboard.py --results {args.output}")


if __name__ == "__main__":
    main()
