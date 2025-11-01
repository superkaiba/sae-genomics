#!/usr/bin/env python3
"""
Main training script for sparse autoencoders on tahoe x1.

Example usage:
    python scripts/train_sae.py --config configs/training/default.yaml
"""

import argparse
import yaml
from pathlib import Path
from rich.console import Console

console = Console()


def load_config(config_path: Path) -> dict:
    """Load configuration from YAML file."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description="Train SAE on tahoe x1 activations")
    parser.add_argument(
        "--config",
        type=Path,
        default="configs/training/default.yaml",
        help="Path to training configuration file",
    )
    parser.add_argument(
        "--model-config",
        type=Path,
        default=None,
        help="Path to model configuration file",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default="./results",
        help="Output directory for results",
    )
    parser.add_argument(
        "--resume",
        type=Path,
        default=None,
        help="Resume training from checkpoint",
    )

    args = parser.parse_args()

    console.print("[bold]SAE Training Script[/bold]")
    console.print(f"Config: {args.config}")

    # Load configuration
    config = load_config(args.config)

    if args.model_config:
        model_config = load_config(args.model_config)
    else:
        model_config_path = config.get("training", {}).get("model_config")
        if model_config_path:
            model_config = load_config(Path(model_config_path))
        else:
            console.print("[red]No model configuration specified[/red]")
            return

    console.print(f"Model: {model_config.get('model', {}).get('name')}")
    console.print(f"Output directory: {args.output_dir}")

    # TODO: Implement training pipeline
    # 1. Load tahoe x1 model
    # 2. Extract activations or set up activation caching
    # 3. Initialize SAE
    # 4. Train SAE using SAELens
    # 5. Save checkpoints and metrics

    console.print("\n[bold red]Training pipeline not yet implemented[/bold red]")
    console.print("\nImplementation steps:")
    console.print("1. Create TahoeModelAdapter in src/sae_genomics/models/")
    console.print("2. Create SAETrainer in src/sae_genomics/training/")
    console.print("3. Set up data loading and activation extraction")
    console.print("4. Integrate with SAELens for training")
    console.print("5. Add experiment tracking (W&B)")


if __name__ == "__main__":
    main()
