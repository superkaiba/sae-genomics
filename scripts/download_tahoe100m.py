#!/usr/bin/env python3
"""
Download the Tahoe-100M dataset from HuggingFace.

This dataset contains 95.6M single-cell samples and is approximately 429 GB.
Uses HuggingFace datasets library to download the expression_data subset.
"""

import argparse
from pathlib import Path
from rich.console import Console
from rich.progress import Progress
from datasets import load_dataset
import subprocess

console = Console()


def main():
    parser = argparse.ArgumentParser(description="Download Tahoe-100M dataset")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/tahoe_100m"),
        help="Output directory for dataset cache",
    )
    parser.add_argument(
        "--subset",
        type=str,
        default="expression_data",
        help="Dataset subset to download (default: expression_data)",
    )
    parser.add_argument(
        "--num-proc",
        type=int,
        default=4,
        help="Number of processes for downloading",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    console.print("[bold green]Downloading Tahoe-100M Dataset[/bold green]")
    console.print(f"Repository: tahoebio/Tahoe-100M")
    console.print(f"Subset: {args.subset}")
    console.print(f"Cache directory: {args.output_dir}")
    console.print(f"Expected size: ~429 GB for expression_data")
    console.print("\n[yellow]Note: This is a large download and may take several hours[/yellow]")
    console.print("[yellow]The dataset will be cached in HuggingFace cache and saved to disk[/yellow]\n")

    try:
        console.print(f"[blue]Loading dataset from HuggingFace...[/blue]")
        console.print(f"[blue]This will download {args.subset} subset with 95.6M rows[/blue]\n")

        # Load the dataset (will download and cache)
        # Use streaming=False to download the full dataset
        dataset = load_dataset(
            "tahoebio/Tahoe-100M",
            args.subset,
            split="train",
            cache_dir=str(args.output_dir / "cache"),
            num_proc=args.num_proc,
        )

        console.print(f"\n[green]✓ Dataset downloaded successfully![/green]")
        console.print(f"[green]Dataset info:[/green]")
        console.print(f"  Rows: {len(dataset):,}")
        console.print(f"  Columns: {dataset.column_names}")
        console.print(f"  Features: {dataset.features}")

        # Save dataset to disk in arrow format for fast loading
        output_path = args.output_dir / f"{args.subset}"
        console.print(f"\n[blue]Saving dataset to {output_path}...[/blue]")
        dataset.save_to_disk(str(output_path))

        console.print(f"[green]✓ Dataset saved to disk[/green]")

        # Show size
        result = subprocess.run(
            ["du", "-sh", str(args.output_dir)],
            capture_output=True,
            text=True,
        )
        console.print(f"\n[bold]Total size:[/bold] {result.stdout.strip()}")

        console.print("\n[bold green]Download complete![/bold green]")
        console.print(f"\nTo load this dataset in Python:")
        console.print(f"  from datasets import load_from_disk")
        console.print(f"  dataset = load_from_disk('{output_path}')")

    except Exception as e:
        console.print(f"[red]✗ Failed to download: {e}[/red]")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
