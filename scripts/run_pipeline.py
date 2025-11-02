#!/usr/bin/env python3
"""
Run the complete SAE genomics pipeline end-to-end.

This master script orchestrates the entire pipeline:
1. Train SAE on Tahoe X1 activations
2. Extract top activating genes for each feature
3. Create visualization dashboard

Usage:
    python scripts/run_pipeline.py \\
        --data path/to/data.h5ad \\
        --model-size 70m \\
        --output results/my_experiment \\
        --steps 10000

For quick testing:
    python scripts/run_pipeline.py \\
        --data path/to/data.h5ad \\
        --model-size 70m \\
        --output results/test \\
        --steps 1000 \\
        --max-cells 1000
"""

import argparse
import subprocess
import sys
from pathlib import Path
from datetime import datetime

from rich.console import Console
from rich.panel import Panel

console = Console()


def run_command(cmd: list[str], description: str) -> bool:
    """Run a command and handle errors.

    Args:
        cmd: Command and arguments
        description: Description of what the command does

    Returns:
        True if successful, False otherwise
    """
    console.print(f"\n[bold cyan]{description}[/bold cyan]")
    console.print(f"[dim]Command: {' '.join(cmd)}[/dim]")

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=False,
            text=True,
        )
        return True
    except subprocess.CalledProcessError as e:
        console.print(f"[bold red]Error: Command failed with exit code {e.returncode}[/bold red]")
        return False
    except KeyboardInterrupt:
        console.print("\n[yellow]Pipeline interrupted by user[/yellow]")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Run complete SAE genomics pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full run with 70M model
  python scripts/run_pipeline.py --data data.h5ad --output results/exp1 --steps 10000

  # Quick test run
  python scripts/run_pipeline.py --data data.h5ad --output results/test --steps 500 --max-cells 500

  # Large model
  python scripts/run_pipeline.py --data data.h5ad --model-size 3b --output results/large --steps 50000
        """,
    )

    # Data arguments
    parser.add_argument(
        "--data",
        type=Path,
        required=True,
        help="Path to AnnData file (.h5ad) with single-cell data",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output directory for all results",
    )

    # Model arguments
    parser.add_argument(
        "--model-size",
        type=str,
        default="70m",
        choices=["70m", "1.3b", "3b"],
        help="Tahoe X1 model size (default: 70m)",
    )
    parser.add_argument(
        "--model-config",
        type=Path,
        default=None,
        help="Path to model config YAML (optional)",
    )

    # Training arguments
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Path to training config YAML (optional)",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=10000,
        help="Training steps (default: 10000)",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=None,
        help="Maximum cells to use (for testing)",
    )

    # Feature extraction arguments
    parser.add_argument(
        "--top-k-genes",
        type=int,
        default=100,
        help="Number of top genes per feature (default: 100)",
    )
    parser.add_argument(
        "--top-k-cells",
        type=int,
        default=100,
        help="Number of top cells per feature (default: 100)",
    )

    # Dashboard arguments
    parser.add_argument(
        "--top-n-features",
        type=int,
        default=100,
        help="Number of features to display in dashboard (default: 100)",
    )

    # Execution control
    parser.add_argument(
        "--device",
        type=str,
        default="auto",
        choices=["auto", "cuda", "cpu"],
        help="Device to use (default: auto)",
    )
    parser.add_argument(
        "--skip-training",
        action="store_true",
        help="Skip training step (use existing checkpoint)",
    )
    parser.add_argument(
        "--skip-extraction",
        action="store_true",
        help="Skip feature extraction step",
    )
    parser.add_argument(
        "--skip-dashboard",
        action="store_true",
        help="Skip dashboard creation step",
    )
    parser.add_argument(
        "--gene-id-key",
        type=str,
        default="ensembl_id",
        help="Key in adata.var for gene IDs",
    )

    args = parser.parse_args()

    # Validate
    if not args.data.exists():
        console.print(f"[red]Error: Data file not found: {args.data}[/red]")
        sys.exit(1)

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Print banner
    start_time = datetime.now()
    banner = f"""
╔══════════════════════════════════════════════════════════╗
║                                                          ║
║            🧬 SAE GENOMICS PIPELINE 🧬                   ║
║                                                          ║
║   Train Sparse Autoencoders on Tahoe X1                ║
║   Extract Interpretable Features                        ║
║   Visualize Biological Insights                         ║
║                                                          ║
╚══════════════════════════════════════════════════════════╝

Configuration:
  Data: {args.data}
  Model: Tahoe X1 {args.model_size}
  Output: {args.output}
  Training Steps: {args.steps:,}
  Max Cells: {args.max_cells or 'All'}
  Device: {args.device}

Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}
"""
    console.print(Panel(banner, border_style="cyan"))

    # Track success of each step
    results = {"training": None, "extraction": None, "dashboard": None}

    # STEP 1: Train SAE
    if not args.skip_training:
        console.print("\n" + "=" * 70)
        console.print("[bold magenta]STEP 1/3: TRAINING SPARSE AUTOENCODER[/bold magenta]")
        console.print("=" * 70)

        cmd = [
            sys.executable,
            "scripts/01_train_sae.py",
            "--data",
            str(args.data),
            "--model-size",
            args.model_size,
            "--output",
            str(args.output),
            "--steps",
            str(args.steps),
            "--device",
            args.device,
            "--gene-id-key",
            args.gene_id_key,
        ]

        if args.max_cells:
            cmd.extend(["--max-cells", str(args.max_cells)])
        if args.model_config:
            cmd.extend(["--model-config", str(args.model_config)])
        if args.config:
            cmd.extend(["--config", str(args.config)])

        results["training"] = run_command(cmd, "Training SAE on Tahoe X1 activations...")

        if not results["training"]:
            console.print("\n[bold red]Training failed. Aborting pipeline.[/bold red]")
            sys.exit(1)
    else:
        console.print("\n[yellow]Skipping training step (--skip-training)[/yellow]")
        results["training"] = "skipped"

    # STEP 2: Extract features
    if not args.skip_extraction:
        console.print("\n" + "=" * 70)
        console.print("[bold magenta]STEP 2/3: EXTRACTING FEATURE ACTIVATIONS[/bold magenta]")
        console.print("=" * 70)

        sae_checkpoint = args.output / "sae_final.pt"
        if not sae_checkpoint.exists():
            console.print(f"[red]Error: SAE checkpoint not found: {sae_checkpoint}[/red]")
            console.print("Please run training first or check output directory")
            sys.exit(1)

        cmd = [
            sys.executable,
            "scripts/02_extract_features.py",
            "--sae",
            str(sae_checkpoint),
            "--data",
            str(args.data),
            "--output",
            str(args.output),
            "--model-size",
            args.model_size,
            "--top-k-genes",
            str(args.top_k_genes),
            "--top-k-cells",
            str(args.top_k_cells),
            "--device",
            args.device,
            "--gene-id-key",
            args.gene_id_key,
            "--use-saved-activations",  # Reuse activations from training
        ]

        if args.max_cells:
            cmd.extend(["--max-cells", str(args.max_cells)])

        results["extraction"] = run_command(
            cmd, "Extracting top activating genes and cells..."
        )

        if not results["extraction"]:
            console.print("\n[bold red]Feature extraction failed. Aborting pipeline.[/bold red]")
            sys.exit(1)
    else:
        console.print("\n[yellow]Skipping feature extraction step (--skip-extraction)[/yellow]")
        results["extraction"] = "skipped"

    # STEP 3: Create dashboard
    if not args.skip_dashboard:
        console.print("\n" + "=" * 70)
        console.print("[bold magenta]STEP 3/3: CREATING VISUALIZATION DASHBOARD[/bold magenta]")
        console.print("=" * 70)

        cmd = [
            sys.executable,
            "scripts/03_create_dashboard.py",
            "--results",
            str(args.output),
            "--top-n",
            str(args.top_n_features),
        ]

        results["dashboard"] = run_command(cmd, "Creating interactive HTML dashboard...")

        if not results["dashboard"]:
            console.print("\n[yellow]Warning: Dashboard creation failed[/yellow]")
    else:
        console.print("\n[yellow]Skipping dashboard creation step (--skip-dashboard)[/yellow]")
        results["dashboard"] = "skipped"

    # Final summary
    end_time = datetime.now()
    duration = end_time - start_time

    summary = f"""
╔══════════════════════════════════════════════════════════╗
║                                                          ║
║              ✅ PIPELINE COMPLETED ✅                     ║
║                                                          ║
╚══════════════════════════════════════════════════════════╝

Results:
  ✓ Training: {results['training']}
  ✓ Feature Extraction: {results['extraction']}
  ✓ Dashboard: {results['dashboard']}

Duration: {duration}
Output Directory: {args.output}

Files Created:
  - {args.output}/sae_final.pt (trained SAE)
  - {args.output}/activations.pt (Tahoe activations)
  - {args.output}/feature_gene_associations.json
  - {args.output}/feature_cell_associations.json
  - {args.output}/dashboard.html

Next Steps:
  1. Open dashboard: open {args.output}/dashboard.html
  2. Analyze specific features
  3. Run biological validation (optional)

For biological validation:
  python scripts/validate_features.py --results {args.output}
"""
    console.print(Panel(summary, border_style="green"))

    # Open dashboard automatically (optional)
    dashboard_path = args.output / "dashboard.html"
    if dashboard_path.exists():
        console.print(f"\n[bold]Dashboard ready at:[/bold] {dashboard_path}")
        console.print("You can open it with: [cyan]open {dashboard_path}[/cyan]")


if __name__ == "__main__":
    main()
