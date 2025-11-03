#!/usr/bin/env python3
"""
Validate SAE features against biological databases.

This script demonstrates the end-to-end validation workflow:
1. Load feature-gene associations from SAE
2. Parse local validation databases
3. Run enrichment analysis
4. Save results

Usage:
    python scripts/validate_sae_features.py --features results/exp1/feature_gene_associations.json
    python scripts/validate_sae_features.py --gene-list IL6,TNF,CXCL8,CCL2,IL1B
"""

import argparse
import json
from pathlib import Path
from typing import List, Dict, Set
from collections import defaultdict

from rich.console import Console
from rich.table import Table
from rich.progress import track

from sae_genomics.validation import LocalValidationClient

console = Console()


def load_feature_genes(feature_file: Path) -> Dict[str, List[str]]:
    """Load feature-gene associations from JSON file.

    Args:
        feature_file: Path to feature_gene_associations.json

    Returns:
        Dict mapping feature IDs to lists of gene symbols
    """
    with open(feature_file) as f:
        data = json.load(f)

    feature_genes = {}
    for feature_id, feature_data in data.items():
        # Extract top gene symbols
        top_genes = feature_data.get('top_genes', [])
        if isinstance(top_genes, list) and len(top_genes) > 0:
            # Handle both formats: [["GENE1", score], ...] or ["GENE1", ...]
            if isinstance(top_genes[0], list):
                genes = [g[0] for g in top_genes]
            else:
                genes = top_genes

            feature_genes[feature_id] = genes

    return feature_genes


def load_client(databases_dir: Path) -> LocalValidationClient:
    """Load local validation client.

    Args:
        databases_dir: Directory containing downloaded databases

    Returns:
        LocalValidationClient instance
    """
    console.print(f"[blue]Loading validation databases from {databases_dir}...[/blue]\n")

    client = LocalValidationClient(
        databases_dir=databases_dir,
        load_on_init=True,
        verbose=True,
    )

    return client


def validate_features(
    feature_genes: Dict[str, List[str]],
    client: LocalValidationClient,
    top_n_features: int = 10,
    max_genes_per_feature: int = 50,
    output_dir: Path = None,
) -> Dict:
    """Run enrichment analysis for features.

    Args:
        feature_genes: Dict mapping feature IDs to gene lists
        client: LocalValidationClient instance
        top_n_features: Number of top features to validate
        max_genes_per_feature: Maximum genes to use per feature
        output_dir: Output directory for results

    Returns:
        Dict with validation results
    """
    # Get database stats
    stats = client.get_database_stats()
    console.print(f"\n[cyan]Background: {stats['total_genes']} genes[/cyan]")

    # Validate each feature
    all_results = {}

    feature_ids = list(feature_genes.keys())[:top_n_features]

    for feature_id in track(feature_ids, description="Validating features"):
        genes = set(feature_genes[feature_id][:max_genes_per_feature])

        if len(genes) < 2:
            continue

        # Run enrichment across all databases
        enrichment_results = client.validate_gene_set(genes)

        # Convert to dict format and limit to top 10 per database
        feature_results = {}
        for db_name, results in enrichment_results.items():
            feature_results[db_name] = [
                r.to_dict() for r in results[:10]
            ]

        all_results[feature_id] = {
            'genes': list(genes),
            'n_genes': len(genes),
            'enrichment': feature_results,
        }

    # Save results
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / 'validation_results.json'

        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=2)

        console.print(f"\n[green]✓ Results saved to {output_file}[/green]")

    return all_results


def print_summary(results: Dict):
    """Print summary table of validation results.

    Args:
        results: Validation results dict
    """
    console.print("\n[bold]Validation Summary[/bold]\n")

    for feature_id, feature_data in results.items():
        enrichment = feature_data['enrichment']

        # Count significant terms across all databases
        db_counts = {
            db_name: len(terms)
            for db_name, terms in enrichment.items()
        }

        console.print(f"[cyan]Feature {feature_id}:[/cyan] {feature_data['n_genes']} genes")

        # Print counts for each database
        for db_name, count in db_counts.items():
            console.print(f"  {db_name}: {count}")

        # Show top result if any
        for db_name, terms in enrichment.items():
            if len(terms) > 0:
                top = terms[0]
                console.print(f"  Top: {top['term_name'][:60]} (p={top['p_value_adjusted']:.2e}, {db_name})")
                break

        console.print()


def main():
    parser = argparse.ArgumentParser(description="Validate SAE features")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--features',
        type=Path,
        help='Path to feature_gene_associations.json from SAE extraction'
    )
    group.add_argument(
        '--gene-list',
        type=str,
        help='Comma-separated list of gene symbols (e.g., "IL6,TNF,CXCL8")'
    )

    parser.add_argument(
        '--databases-dir',
        type=Path,
        default=Path('data/validation_databases'),
        help='Directory containing downloaded validation databases'
    )

    parser.add_argument(
        '--output-dir',
        type=Path,
        help='Output directory for validation results'
    )

    parser.add_argument(
        '--top-n',
        type=int,
        default=10,
        help='Number of top features to validate (default: 10)'
    )

    args = parser.parse_args()

    # Load features
    if args.features:
        console.print(f"[blue]Loading features from {args.features}...[/blue]")
        feature_genes = load_feature_genes(args.features)
        console.print(f"[green]✓ Loaded {len(feature_genes)} features[/green]\n")
    else:
        # Single gene list
        genes = [g.strip() for g in args.gene_list.split(',')]
        feature_genes = {'query': genes}
        console.print(f"[green]✓ Testing {len(genes)} genes[/green]\n")

    # Load validation client
    try:
        client = load_client(args.databases_dir)
    except Exception as e:
        console.print(f"[red]✗ Failed to load databases: {e}[/red]")
        console.print("[yellow]Download databases first:[/yellow]")
        console.print("  python scripts/download_validation_databases.py --tier 1")
        return 1

    # Check if any databases loaded
    stats = client.get_database_stats()
    if stats['total_genes'] == 0:
        console.print("[red]✗ No databases found. Download databases first:[/red]")
        console.print("  python scripts/download_validation_databases.py --tier 1")
        return 1

    # Run validation
    console.print(f"\n[bold]Running enrichment analysis...[/bold]\n")
    results = validate_features(
        feature_genes=feature_genes,
        client=client,
        top_n_features=args.top_n,
        output_dir=args.output_dir,
    )

    # Print summary
    print_summary(results)

    console.print("[green]✓ Validation complete![/green]")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
