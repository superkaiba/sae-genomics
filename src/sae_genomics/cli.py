"""Command-line interface for SAE Genomics."""

import typer
from pathlib import Path
from typing import Optional
from rich.console import Console
from rich.progress import track

app = typer.Typer(help="SAE Genomics: Train and validate sparse autoencoders on tahoe x1")
console = Console()


@app.command()
def train(
    config: Path = typer.Option(
        "configs/training/default.yaml",
        "--config",
        "-c",
        help="Path to training configuration file",
    ),
    model_config: Optional[Path] = typer.Option(
        None,
        "--model-config",
        "-m",
        help="Path to model configuration file (overrides config)",
    ),
    output_dir: Path = typer.Option(
        "./results",
        "--output",
        "-o",
        help="Output directory for results",
    ),
    resume: Optional[Path] = typer.Option(
        None,
        "--resume",
        "-r",
        help="Resume training from checkpoint",
    ),
):
    """Train a sparse autoencoder on tahoe x1 activations."""
    console.print(f"[bold green]Starting SAE training[/bold green]")
    console.print(f"Config: {config}")
    console.print(f"Output directory: {output_dir}")

    # TODO: Implement training pipeline
    # from sae_genomics.training import SAETrainer
    # trainer = SAETrainer(config, model_config, output_dir)
    # trainer.train()

    console.print("[bold red]Training pipeline not yet implemented[/bold red]")
    console.print("Please implement the training pipeline in src/sae_genomics/training/trainer.py")


@app.command()
def validate(
    sae_path: Path = typer.Argument(..., help="Path to trained SAE checkpoint"),
    config: Path = typer.Option(
        "configs/validation/databases.yaml",
        "--config",
        "-c",
        help="Path to validation configuration file",
    ),
    output_dir: Path = typer.Option(
        "./results/validation",
        "--output",
        "-o",
        help="Output directory for validation results",
    ),
    feature_ids: Optional[str] = typer.Option(
        None,
        "--features",
        "-f",
        help="Comma-separated feature IDs to validate (default: all)",
    ),
):
    """Run biological validation on trained SAE features."""
    console.print(f"[bold green]Starting biological validation[/bold green]")
    console.print(f"SAE checkpoint: {sae_path}")
    console.print(f"Config: {config}")
    console.print(f"Output directory: {output_dir}")

    # TODO: Implement validation pipeline
    # from sae_genomics.validation import FeatureValidator
    # validator = FeatureValidator(config)
    # results = validator.validate(sae_path, feature_ids)
    # validator.save_results(results, output_dir)

    console.print("[bold red]Validation pipeline not yet implemented[/bold red]")
    console.print("Please implement the validation pipeline in src/sae_genomics/validation/")


@app.command()
def export_results(
    results_dir: Path = typer.Argument(..., help="Directory containing training results"),
    output_format: str = typer.Option(
        "all",
        "--format",
        "-f",
        help="Export format: checkpoint, features, metrics, all",
    ),
    output_path: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output path (default: results_dir/export)",
    ),
):
    """Export trained SAE and results to various formats."""
    console.print(f"[bold green]Exporting results[/bold green]")
    console.print(f"Results directory: {results_dir}")
    console.print(f"Format: {output_format}")

    if output_path is None:
        output_path = results_dir / "export"

    console.print(f"Output path: {output_path}")

    # TODO: Implement export functionality
    console.print("[bold red]Export functionality not yet implemented[/bold red]")


if __name__ == "__main__":
    app()
