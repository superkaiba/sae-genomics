#!/usr/bin/env python3
"""
Convert feature_gene_associations.json from Ensembl IDs to gene symbols.

This script fixes the mismatch between extracted features (which use Ensembl IDs)
and validation databases (which use gene symbols).

Usage:
    python scripts/convert_ensembl_to_symbols.py results/quick_test
"""

import argparse
import json
from pathlib import Path
import sys

import scanpy as sc
from rich.console import Console
from rich.progress import track

console = Console()


def load_ensembl_to_symbol_mapping(adata_path: Path) -> dict:
    """Load Ensembl ID to gene symbol mapping from AnnData.

    Args:
        adata_path: Path to .h5ad file

    Returns:
        Dict mapping Ensembl IDs to gene symbols
    """
    console.print(f"[blue]Loading gene mapping from {adata_path}...[/blue]")

    adata = sc.read_h5ad(adata_path)

    # Create mapping: Ensembl ID → gene symbol
    ensembl_to_symbol = {}

    for gene_symbol in adata.var.index:
        ensembl_id = adata.var.loc[gene_symbol, 'ensembl_id']
        ensembl_to_symbol[ensembl_id] = gene_symbol

    console.print(f"[green]✓ Loaded {len(ensembl_to_symbol)} gene mappings[/green]")

    return ensembl_to_symbol


def convert_feature_associations(
    input_file: Path,
    output_file: Path,
    ensembl_to_symbol: dict,
):
    """Convert feature associations from Ensembl IDs to gene symbols.

    Args:
        input_file: Input JSON file with Ensembl IDs
        output_file: Output JSON file with gene symbols
        ensembl_to_symbol: Mapping dict
    """
    console.print(f"\n[blue]Loading feature associations from {input_file}...[/blue]")

    with open(input_file) as f:
        data = json.load(f)

    console.print(f"[green]✓ Loaded {len(data)} features[/green]")

    # Convert each feature's genes
    n_converted = 0
    n_unmapped = 0
    n_total_genes = 0

    converted_data = {}

    console.print("\n[blue]Converting Ensembl IDs to gene symbols...[/blue]")

    for feature_id, feature_data in track(data.items(), description="Converting"):
        top_genes = feature_data.get('top_genes', [])

        # Convert gene IDs to symbols
        converted_top_genes = []
        for gene_entry in top_genes:
            if isinstance(gene_entry, list):
                ensembl_id, score = gene_entry
            else:
                ensembl_id = gene_entry
                score = 0.0

            n_total_genes += 1

            # Convert to symbol
            if ensembl_id in ensembl_to_symbol:
                gene_symbol = ensembl_to_symbol[ensembl_id]
                converted_top_genes.append([gene_symbol, score])
                n_converted += 1
            elif ensembl_id.startswith('ENSG'):
                # Unmapped Ensembl ID
                n_unmapped += 1
            elif ensembl_id in ['<cls>', '<pad>', '<eos>']:
                # Special tokens - skip
                pass
            else:
                # Already a symbol or unknown
                converted_top_genes.append([ensembl_id, score])

        # Convert gene_scores dict
        gene_scores = feature_data.get('gene_scores', {})
        converted_gene_scores = {}
        for ensembl_id, score in gene_scores.items():
            if ensembl_id in ensembl_to_symbol:
                gene_symbol = ensembl_to_symbol[ensembl_id]
                converted_gene_scores[gene_symbol] = score
            elif not ensembl_id.startswith('ENSG'):
                # Already a symbol
                converted_gene_scores[ensembl_id] = score

        # Store converted data
        converted_data[feature_id] = {
            'top_genes': converted_top_genes,
            'gene_scores': converted_gene_scores,
            'n_activations': feature_data.get('n_activations', 0),
        }

    # Save converted data
    console.print(f"\n[blue]Saving converted associations to {output_file}...[/blue]")
    with open(output_file, 'w') as f:
        json.dump(converted_data, f, indent=2)

    console.print(f"[green]✓ Saved converted associations[/green]")
    console.print(f"\nConversion stats:")
    console.print(f"  Total genes: {n_total_genes}")
    console.print(f"  Converted: {n_converted}")
    console.print(f"  Unmapped: {n_unmapped}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert feature associations from Ensembl IDs to gene symbols"
    )
    parser.add_argument(
        'results_dir',
        type=Path,
        help='Results directory containing feature_gene_associations.json'
    )
    parser.add_argument(
        '--data',
        type=Path,
        default=Path('data/pbmc3k_test.h5ad'),
        help='Path to original .h5ad file for gene mapping (default: data/pbmc3k_test.h5ad)'
    )
    parser.add_argument(
        '--output',
        type=Path,
        help='Output file path (default: <results_dir>/feature_gene_associations_symbols.json)'
    )

    args = parser.parse_args()

    # Find input file
    input_file = args.results_dir / 'feature_gene_associations.json'
    if not input_file.exists():
        console.print(f"[red]✗ File not found: {input_file}[/red]")
        return 1

    # Set output file
    if args.output is None:
        args.output = args.results_dir / 'feature_gene_associations_symbols.json'

    console.print("[bold]Feature Association Converter[/bold]")
    console.print(f"Input: {input_file}")
    console.print(f"Data: {args.data}")
    console.print(f"Output: {args.output}")

    # Load mapping
    try:
        ensembl_to_symbol = load_ensembl_to_symbol_mapping(args.data)
    except Exception as e:
        console.print(f"[red]✗ Failed to load gene mapping: {e}[/red]")
        return 1

    # Convert
    try:
        convert_feature_associations(
            input_file=input_file,
            output_file=args.output,
            ensembl_to_symbol=ensembl_to_symbol,
        )
    except Exception as e:
        console.print(f"[red]✗ Conversion failed: {e}[/red]")
        import traceback
        traceback.print_exc()
        return 1

    console.print("\n[green]✓ Conversion complete![/green]")
    console.print(f"\n[yellow]Next step: Re-run validation with converted features:[/yellow]")
    console.print(f"  python scripts/validate_sae_features.py --features {args.output}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
