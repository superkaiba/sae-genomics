#!/usr/bin/env python3
"""
Step 4: Validate SAE features against biological databases.

This script takes feature-gene associations from step 02 and validates
them against multiple biological databases using statistical enrichment
analysis.

Usage:
    python scripts/04_validate_features.py \\
        --features results/exp1/feature_gene_associations.json \\
        --output results/exp1/validation

Pipeline integration:
    Step 1: 01_train_sae.py → Train SAE
    Step 2: 02_extract_features.py → Extract feature-gene associations
    Step 3: 03_create_dashboard.py → Create visualization
    Step 4: 04_validate_features.py → Validate features (THIS SCRIPT)
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List

from rich.console import Console
from rich.table import Table
from rich.progress import track

from sae_genomics.validation import LocalValidationClient

console = Console()


def load_feature_associations(features_file: Path) -> Dict[str, List[str]]:
    """Load feature-gene associations from JSON file.

    Args:
        features_file: Path to feature_gene_associations.json from step 02

    Returns:
        Dict mapping feature IDs to lists of gene symbols
    """
    console.print(f"[blue]Loading features from {features_file}...[/blue]")

    with open(features_file) as f:
        data = json.load(f)

    feature_genes = {}

    for feature_id, feature_data in data.items():
        # Extract gene symbols from top genes
        # Format: {"top_genes": [["GENE1", score], ["GENE2", score], ...]}
        top_genes = feature_data.get('top_genes', [])

        if not top_genes:
            continue

        # Handle both formats: [[gene, score], ...] or [gene, ...]
        if isinstance(top_genes[0], list):
            genes = [g[0] for g in top_genes]
        else:
            genes = top_genes

        feature_genes[feature_id] = genes

    console.print(f"[green]✓ Loaded {len(feature_genes)} features[/green]\n")

    return feature_genes


def validate_features(
    feature_genes: Dict[str, List[str]],
    client: LocalValidationClient,
    max_features: int = None,
    max_genes_per_feature: int = 50,
) -> Dict:
    """Run enrichment analysis for all features.

    Args:
        feature_genes: Dict mapping feature IDs to gene lists
        client: LocalValidationClient instance
        max_features: Maximum number of features to validate (default: all)
        max_genes_per_feature: Max genes per feature for enrichment

    Returns:
        Dict with validation results
    """
    # Determine features to validate
    feature_ids = list(feature_genes.keys())

    if max_features:
        feature_ids = feature_ids[:max_features]

    console.print(f"[cyan]Validating {len(feature_ids)} features...[/cyan]\n")

    # Get database stats
    stats = client.get_database_stats()
    console.print(f"[cyan]Background: {stats['total_genes']} genes across {len(stats['databases'])} databases[/cyan]\n")

    # Validate each feature
    all_results = {}
    n_with_enrichments = 0

    for feature_id in track(feature_ids, description="Running enrichment analysis"):
        # Get top genes for this feature
        genes = set(feature_genes[feature_id][:max_genes_per_feature])

        if len(genes) < 2:
            continue

        # Run enrichment across all databases
        enrichment_results = client.validate_gene_set(genes)

        # Convert to serializable format and count enrichments
        feature_enrichments = {}
        has_enrichments = False

        for db_name, results in enrichment_results.items():
            if results:
                has_enrichments = True
                feature_enrichments[db_name] = [r.to_dict() for r in results[:20]]  # Top 20 per database

        if has_enrichments:
            n_with_enrichments += 1

        all_results[feature_id] = {
            'genes': list(genes),
            'n_genes': len(genes),
            'enrichments': feature_enrichments,
            'n_enrichments_total': sum(len(r) for r in enrichment_results.values()),
        }

    console.print(f"\n[green]✓ {n_with_enrichments}/{len(feature_ids)} features have significant enrichments[/green]\n")

    return all_results


def print_summary(results: Dict, top_n: int = 10):
    """Print summary of validation results.

    Args:
        results: Validation results dict
        top_n: Number of top features to display
    """
    console.print("[bold]Validation Summary[/bold]\n")

    # Sort features by total enrichments
    sorted_features = sorted(
        results.items(),
        key=lambda x: x[1]['n_enrichments_total'],
        reverse=True
    )

    # Create summary table
    table = Table(title=f"Top {top_n} Features by Enrichments")
    table.add_column("Feature", style="cyan")
    table.add_column("Genes", justify="right")
    table.add_column("Total Enrichments", justify="right")
    table.add_column("Top Result", style="green")

    for feature_id, data in sorted_features[:top_n]:
        n_genes = data['n_genes']
        n_enrichments = data['n_enrichments_total']

        # Get top enrichment result
        top_result = ""
        enrichments = data['enrichments']

        for db_name, results in enrichments.items():
            if results:
                top = results[0]
                term_name = top['term_name'][:40]
                p_value = top['p_value_adjusted']
                top_result = f"{term_name} (p={p_value:.2e}, {db_name})"
                break

        table.add_row(
            feature_id,
            str(n_genes),
            str(n_enrichments),
            top_result
        )

    console.print(table)
    console.print()


def save_results(results: Dict, output_dir: Path):
    """Save validation results to files.

    Args:
        results: Validation results dict
        output_dir: Output directory
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save full results as JSON
    results_file = output_dir / 'feature_validations.json'
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)

    console.print(f"[green]✓ Full results saved to {results_file}[/green]")

    # Create summary CSV for easy viewing
    import csv

    summary_file = output_dir / 'validation_summary.csv'
    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'feature_id',
            'n_genes',
            'n_enrichments_total',
            'n_gencc',
            'n_hpo_phenotypes',
            'n_hpo_diseases',
            'n_gwas',
            'n_orphanet',
            'top_term',
            'top_p_value',
            'top_database',
        ])

        for feature_id, data in results.items():
            enrichments = data['enrichments']

            # Count per database
            counts = {
                'gencc': len(enrichments.get('gencc', [])),
                'hpo_phenotypes': len(enrichments.get('hpo_phenotypes', [])),
                'hpo_diseases': len(enrichments.get('hpo_diseases', [])),
                'gwas': len(enrichments.get('gwas', [])),
                'orphanet': len(enrichments.get('orphanet', [])),
            }

            # Get top result
            top_term = ''
            top_p = ''
            top_db = ''

            for db_name, results in enrichments.items():
                if results:
                    top_term = results[0]['term_name']
                    top_p = results[0]['p_value_adjusted']
                    top_db = db_name
                    break

            writer.writerow([
                feature_id,
                data['n_genes'],
                data['n_enrichments_total'],
                counts['gencc'],
                counts['hpo_phenotypes'],
                counts['hpo_diseases'],
                counts['gwas'],
                counts['orphanet'],
                top_term,
                top_p,
                top_db,
            ])

    console.print(f"[green]✓ Summary saved to {summary_file}[/green]\n")


def main():
    parser = argparse.ArgumentParser(
        description="Validate SAE features against biological databases",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Validate all features
    python scripts/04_validate_features.py \\
        --features results/exp1/feature_gene_associations.json \\
        --output results/exp1/validation

    # Validate top 100 features only
    python scripts/04_validate_features.py \\
        --features results/exp1/feature_gene_associations.json \\
        --output results/exp1/validation \\
        --max-features 100

    # Use custom database directory
    python scripts/04_validate_features.py \\
        --features results/exp1/feature_gene_associations.json \\
        --output results/exp1/validation \\
        --databases-dir /path/to/databases
        """
    )

    parser.add_argument(
        '--features',
        type=Path,
        required=True,
        help='Path to feature_gene_associations.json from step 02'
    )

    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output directory for validation results'
    )

    parser.add_argument(
        '--databases-dir',
        type=Path,
        default=Path('data/validation_databases'),
        help='Directory containing validation databases (default: data/validation_databases)'
    )

    parser.add_argument(
        '--max-features',
        type=int,
        help='Maximum number of features to validate (default: all)'
    )

    parser.add_argument(
        '--max-genes',
        type=int,
        default=50,
        help='Maximum genes per feature for enrichment (default: 50)'
    )

    parser.add_argument(
        '--summary-top-n',
        type=int,
        default=20,
        help='Number of features to show in summary (default: 20)'
    )

    args = parser.parse_args()

    # Check inputs
    if not args.features.exists():
        console.print(f"[red]✗ Features file not found: {args.features}[/red]")
        return 1

    if not args.databases_dir.exists():
        console.print(f"[red]✗ Databases directory not found: {args.databases_dir}[/red]")
        console.print("[yellow]Download databases first:[/yellow]")
        console.print("  python scripts/download_validation_databases.py --tier 1 2")
        return 1

    # Load feature associations
    feature_genes = load_feature_associations(args.features)

    if not feature_genes:
        console.print("[red]✗ No features found in input file[/red]")
        return 1

    # Load validation client
    console.print("[blue]Loading validation databases...[/blue]\n")

    try:
        client = LocalValidationClient(
            databases_dir=args.databases_dir,
            load_on_init=True,
            verbose=True,
        )
    except Exception as e:
        console.print(f"[red]✗ Failed to load databases: {e}[/red]")
        return 1

    # Check databases loaded
    stats = client.get_database_stats()
    if stats['total_genes'] == 0:
        console.print("[red]✗ No databases loaded[/red]")
        console.print("[yellow]Download databases first:[/yellow]")
        console.print("  python scripts/download_validation_databases.py --tier 1 2")
        return 1

    # Run validation
    console.print(f"\n[bold]Running enrichment analysis...[/bold]\n")

    results = validate_features(
        feature_genes=feature_genes,
        client=client,
        max_features=args.max_features,
        max_genes_per_feature=args.max_genes,
    )

    # Print summary
    print_summary(results, top_n=args.summary_top_n)

    # Save results
    save_results(results, args.output)

    console.print("[bold green]✓ Validation complete![/bold green]")

    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
