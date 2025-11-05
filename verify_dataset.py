#!/usr/bin/env python3
"""Quick verification that Tahoe-100M dataset loads correctly."""

from datasets import load_from_disk
from rich.console import Console
from pathlib import Path

console = Console()

dataset_path = Path("data/tahoe_100m/expression_data")

console.print(f"[blue]Loading dataset from {dataset_path}...[/blue]")

try:
    dataset = load_from_disk(str(dataset_path))

    console.print(f"[green]✓ Dataset loaded successfully![/green]")
    console.print(f"\n[bold]Dataset Info:[/bold]")
    console.print(f"  Total examples: {len(dataset):,}")
    console.print(f"  Columns: {dataset.column_names}")
    console.print(f"  Features: {dataset.features}")

    # Test loading first example
    console.print(f"\n[blue]Testing first example...[/blue]")
    first = dataset[0]
    console.print(f"[green]✓ First example loaded[/green]")
    console.print(f"  Keys: {list(first.keys())}")
    console.print(f"  Genes array length: {len(first['genes'])}")
    console.print(f"  Expressions array length: {len(first['expressions'])}")

    console.print(f"\n[bold green]Dataset verification PASSED! ✓[/bold green]")
    console.print(f"[green]The dataset is ready for SAE training.[/green]")

except Exception as e:
    console.print(f"[red]✗ Error loading dataset: {e}[/red]")
    import traceback
    traceback.print_exc()
