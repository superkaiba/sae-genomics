#!/usr/bin/env python3
"""
Create an interactive dashboard visualizing SAE features.

Usage:
    python scripts/03_create_dashboard.py \\
        --results results/exp1 \\
        --output results/exp1/dashboard.html

This creates an HTML dashboard showing:
- All SAE features
- Top activating genes per feature
- Biological annotations (GO terms, pathways)
- Cell type associations
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List

from rich.console import Console

console = Console()


def create_html_dashboard(
    gene_results: Dict,
    cell_results: Dict,
    summary: Dict,
    output_path: Path,
    top_n_features: int = 50,
):
    """Create an HTML dashboard for SAE features.

    Args:
        gene_results: Feature-gene associations
        cell_results: Feature-cell associations
        summary: Extraction summary
        output_path: Path to save HTML file
        top_n_features: Number of top features to display
    """
    # Sort features by number of activations
    sorted_features = sorted(
        gene_results.items(),
        key=lambda x: x[1].get("n_activations", 0),
        reverse=True,
    )[:top_n_features]

    html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SAE Feature Dashboard</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            background: #0a0a0a;
            color: #e0e0e0;
            line-height: 1.6;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}
        header {{
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
            padding: 40px 20px;
            border-radius: 12px;
            margin-bottom: 30px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.3);
        }}
        h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }}
        .summary {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .summary-card {{
            background: #1a1a2e;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .summary-card h3 {{
            color: #667eea;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 1px;
            margin-bottom: 10px;
        }}
        .summary-card .value {{
            font-size: 2em;
            font-weight: bold;
        }}
        .feature-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .feature-card {{
            background: #1a1a2e;
            border-radius: 12px;
            padding: 20px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.2);
            transition: transform 0.2s, box-shadow 0.2s;
        }}
        .feature-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 8px 24px rgba(102, 126, 234, 0.3);
        }}
        .feature-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 15px;
            padding-bottom: 15px;
            border-bottom: 2px solid #2a2a3e;
        }}
        .feature-id {{
            font-size: 1.5em;
            font-weight: bold;
            color: #667eea;
        }}
        .feature-stats {{
            text-align: right;
            font-size: 0.9em;
            color: #888;
        }}
        .gene-list {{
            margin-bottom: 15px;
        }}
        .gene-list h4 {{
            color: #667eea;
            margin-bottom: 10px;
            font-size: 1em;
        }}
        .gene-item {{
            display: flex;
            justify-content: space-between;
            padding: 8px 12px;
            margin-bottom: 4px;
            background: #0f0f1e;
            border-radius: 6px;
            font-size: 0.9em;
        }}
        .gene-name {{
            color: #e0e0e0;
            font-weight: 500;
        }}
        .gene-score {{
            color: #888;
            font-family: 'Courier New', monospace;
        }}
        .cell-types {{
            margin-top: 15px;
            padding-top: 15px;
            border-top: 1px solid #2a2a3e;
        }}
        .cell-types h4 {{
            color: #764ba2;
            margin-bottom: 10px;
            font-size: 0.9em;
        }}
        .cell-type-badge {{
            display: inline-block;
            background: #16213e;
            color: #e0e0e0;
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 0.85em;
            margin-right: 8px;
            margin-bottom: 8px;
        }}
        .search-box {{
            position: sticky;
            top: 20px;
            z-index: 100;
            margin-bottom: 30px;
        }}
        .search-box input {{
            width: 100%;
            padding: 15px 20px;
            font-size: 1.1em;
            background: #1a1a2e;
            border: 2px solid #2a2a3e;
            border-radius: 8px;
            color: #e0e0e0;
            transition: border-color 0.3s;
        }}
        .search-box input:focus {{
            outline: none;
            border-color: #667eea;
        }}
        .no-results {{
            text-align: center;
            padding: 60px 20px;
            color: #888;
            font-size: 1.2em;
        }}
        footer {{
            text-align: center;
            padding: 40px 20px;
            color: #666;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🧬 SAE Feature Dashboard</h1>
            <p style="color: #888; margin-top: 10px;">Interpretable features from sparse autoencoders on single-cell data</p>
        </header>

        <div class="summary">
            <div class="summary-card">
                <h3>Total Features</h3>
                <div class="value">{summary.get('n_features', 'N/A')}</div>
            </div>
            <div class="summary-card">
                <h3>Active Features</h3>
                <div class="value">{summary.get('active_features', 'N/A')}</div>
            </div>
            <div class="summary-card">
                <h3>Cells Analyzed</h3>
                <div class="value">{summary.get('n_cells_analyzed', 'N/A'):,}</div>
            </div>
            <div class="summary-card">
                <h3>Top Genes/Feature</h3>
                <div class="value">{summary.get('top_k_genes', 'N/A')}</div>
            </div>
        </div>

        <div class="search-box">
            <input type="text" id="searchInput" placeholder="🔍 Search features by gene name or ID...">
        </div>

        <div class="feature-grid" id="featureGrid">
"""

    # Add feature cards
    for feat_id, feat_data in sorted_features:
        top_genes = feat_data.get("top_genes", [])[:10]
        n_activations = feat_data.get("n_activations", 0)

        # Get cell type info if available
        cell_data = cell_results.get(str(feat_id), {})
        cell_types = cell_data.get("cell_types", {})

        html += f"""
            <div class="feature-card" data-feature-id="{feat_id}">
                <div class="feature-header">
                    <div class="feature-id">Feature {feat_id}</div>
                    <div class="feature-stats">
                        {n_activations:,} activations
                    </div>
                </div>

                <div class="gene-list">
                    <h4>Top Activating Genes</h4>
"""

        for gene_name, score in top_genes:
            html += f"""
                    <div class="gene-item">
                        <span class="gene-name">{gene_name}</span>
                        <span class="gene-score">{score:.4f}</span>
                    </div>
"""

        html += """
                </div>
"""

        # Add cell type information if available
        if cell_types:
            sorted_cell_types = sorted(cell_types.items(), key=lambda x: x[1], reverse=True)[
                :5
            ]
            html += """
                <div class="cell-types">
                    <h4>Enriched Cell Types</h4>
"""
            for cell_type, count in sorted_cell_types:
                html += f"""
                    <span class="cell-type-badge">{cell_type} ({count})</span>
"""
            html += """
                </div>
"""

        html += """
            </div>
"""

    html += """
        </div>

        <div class="no-results" id="noResults" style="display: none;">
            No features found matching your search.
        </div>

        <footer>
            Generated with SAE Genomics | Sparse Autoencoders on Tahoe X1
        </footer>
    </div>

    <script>
        // Search functionality
        const searchInput = document.getElementById('searchInput');
        const featureGrid = document.getElementById('featureGrid');
        const noResults = document.getElementById('noResults');
        const featureCards = document.querySelectorAll('.feature-card');

        searchInput.addEventListener('input', (e) => {
            const searchTerm = e.target.value.toLowerCase();
            let visibleCount = 0;

            featureCards.forEach(card => {
                const text = card.textContent.toLowerCase();
                if (text.includes(searchTerm)) {
                    card.style.display = 'block';
                    visibleCount++;
                } else {
                    card.style.display = 'none';
                }
            });

            // Show/hide no results message
            if (visibleCount === 0) {
                featureGrid.style.display = 'none';
                noResults.style.display = 'block';
            } else {
                featureGrid.style.display = 'grid';
                noResults.style.display = 'none';
            }
        });

        // Add smooth scrolling
        document.querySelectorAll('a[href^="#"]').forEach(anchor => {
            anchor.addEventListener('click', function (e) {
                e.preventDefault();
                const target = document.querySelector(this.getAttribute('href'));
                if (target) {
                    target.scrollIntoView({ behavior: 'smooth' });
                }
            });
        });
    </script>
</body>
</html>
"""

    # Write to file
    with open(output_path, "w") as f:
        f.write(html)

    console.print(f"✓ Created HTML dashboard: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Create SAE feature dashboard")
    parser.add_argument(
        "--results",
        type=Path,
        required=True,
        help="Results directory containing feature extractions",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output HTML file path (default: results/dashboard.html)",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=100,
        help="Number of top features to display",
    )

    args = parser.parse_args()

    # Set output path
    if args.output is None:
        args.output = args.results / "dashboard.html"

    console.print("[bold green]Dashboard Creation Pipeline[/bold green]")
    console.print(f"Results directory: {args.results}")
    console.print(f"Output: {args.output}")

    # Step 1: Load results
    console.print("\n[bold]Step 1: Loading extraction results[/bold]")

    gene_results_path = args.results / "feature_gene_associations.json"
    cell_results_path = args.results / "feature_cell_associations.json"
    summary_path = args.results / "extraction_summary.json"

    # Check files exist
    if not gene_results_path.exists():
        console.print(f"[red]Error: {gene_results_path} not found[/red]")
        console.print("Please run 02_extract_features.py first")
        return

    # Load gene results
    with open(gene_results_path) as f:
        gene_results = json.load(f)
    console.print(f"✓ Loaded gene associations for {len(gene_results)} features")

    # Load cell results
    cell_results = {}
    if cell_results_path.exists():
        with open(cell_results_path) as f:
            cell_results = json.load(f)
        console.print(f"✓ Loaded cell associations")
    else:
        console.print("[yellow]  Cell associations not found, skipping[/yellow]")

    # Load summary
    summary = {}
    if summary_path.exists():
        with open(summary_path) as f:
            summary = json.load(f)
        console.print(f"✓ Loaded summary")
    else:
        console.print("[yellow]  Summary not found, using defaults[/yellow]")

    # Step 2: Create dashboard
    console.print("\n[bold]Step 2: Creating HTML dashboard[/bold]")

    create_html_dashboard(
        gene_results=gene_results,
        cell_results=cell_results,
        summary=summary,
        output_path=args.output,
        top_n_features=args.top_n,
    )

    # Final summary
    console.print("\n[bold green]Dashboard Created Successfully![/bold green]")
    console.print(f"Dashboard saved to: {args.output}")
    console.print(f"\nTo view:")
    console.print(f"  open {args.output}")
    console.print(f"\nFeatures:")
    console.print(f"  - Search functionality to filter features by gene")
    console.print(f"  - Top {args.top_n} most active features displayed")
    console.print(f"  - Gene and cell type associations shown")


if __name__ == "__main__":
    main()
