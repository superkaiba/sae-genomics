#!/usr/bin/env python3
"""
Extract top activating genes for each SAE feature.

Usage:
    python scripts/02_extract_features.py \\
        --sae results/exp1/sae_final.pt \\
        --data path/to/data.h5ad \\
        --output results/exp1

This script runs cells through the trained SAE and identifies which genes
most strongly activate each feature.
"""

import argparse
import json
from pathlib import Path

import scanpy as sc
import torch
from rich.console import Console

from sae_genomics.models.model_config import ModelConfig
from sae_genomics.models.tahoe_adapter import TahoeModelAdapter
from sae_genomics.training.trainer import SAETrainer
from sae_genomics.utils.feature_analysis import FeatureAnalyzer

console = Console()

# Configuration constants
EXTRACTION_BATCH_SIZE = 32  # Batch size for activation extraction
EXTRACTION_MAX_LENGTH = 1024  # Maximum sequence length
EXTRACTION_NUM_WORKERS = 4  # Number of dataloader workers


def main():
    parser = argparse.ArgumentParser(description="Extract SAE feature activations")
    parser.add_argument(
        "--sae",
        type=Path,
        required=True,
        help="Path to trained SAE checkpoint (.pt)",
    )
    parser.add_argument(
        "--data",
        type=Path,
        required=True,
        help="Path to AnnData file (.h5ad)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output directory (defaults to SAE directory)",
    )
    parser.add_argument(
        "--model-size",
        type=str,
        default="70m",
        choices=["70m", "1.3b", "3b"],
        help="Tahoe X1 model size",
    )
    parser.add_argument(
        "--top-k-genes",
        type=int,
        default=100,
        help="Number of top genes to extract per feature",
    )
    parser.add_argument(
        "--top-k-cells",
        type=int,
        default=100,
        help="Number of top cells to extract per feature",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=None,
        help="Maximum number of cells to analyze",
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
    parser.add_argument(
        "--use-saved-activations",
        action="store_true",
        help="Use pre-saved activations if available",
    )

    args = parser.parse_args()

    # Set output directory
    if args.output is None:
        args.output = args.sae.parent
    args.output.mkdir(parents=True, exist_ok=True)

    console.print("[bold green]Feature Extraction Pipeline[/bold green]")
    console.print(f"SAE: {args.sae}")
    console.print(f"Data: {args.data}")
    console.print(f"Output: {args.output}")

    # Step 1: Load trained SAE
    console.print("\n[bold]Step 1: Loading trained SAE[/bold]")
    trainer = SAETrainer.load_checkpoint(args.sae, device=args.device)
    console.print(f"✓ Loaded SAE")
    console.print(f"  - d_in: {trainer.d_in}")
    console.print(f"  - d_sae: {trainer.d_sae} features")
    console.print(f"  - Activation: {trainer.activation}")

    # Step 2: Load Tahoe X1 model
    console.print("\n[bold]Step 2: Loading Tahoe X1 model[/bold]")
    tahoe_adapter = TahoeModelAdapter.from_hf(
        model_size=args.model_size,
        device=args.device,
    )
    console.print(f"✓ Loaded {args.model_size} model")

    # Create feature analyzer
    # Get vocab dict from GeneVocab object
    vocab_dict = tahoe_adapter.vocab.get_stoi()
    analyzer = FeatureAnalyzer(
        vocab=vocab_dict,
        reverse_vocab={v: k for k, v in vocab_dict.items()},
    )

    # Step 3: Load or extract activations
    console.print("\n[bold]Step 3: Preparing activations[/bold]")

    activations_path = args.output / "activations.pt"
    cell_metadata = None  # Initialize cell metadata

    if args.use_saved_activations and activations_path.exists():
        console.print(f"Loading saved activations from {activations_path}")
        saved_data = torch.load(activations_path)
        activations = saved_data["activations"]
        gene_ids = saved_data["gene_ids"]
        console.print(f"✓ Loaded activations: {activations.shape}")

        # Load data for cell metadata
        try:
            console.print(f"Loading cell metadata from {args.data}")
            adata = sc.read_h5ad(args.data)
            cell_metadata = adata.obs
            console.print(f"  Loaded cell metadata with {len(cell_metadata.columns)} columns")
        except Exception as e:
            console.print(f"  [yellow]Warning: Could not load cell metadata: {e}[/yellow]")
    else:
        # Load data
        console.print(f"Loading data from {args.data}")
        adata = sc.read_h5ad(args.data)
        console.print(f"✓ Loaded {adata.n_obs} cells × {adata.n_vars} genes")

        # Extract cell metadata before subsampling
        cell_metadata = adata.obs

        # Subset if requested
        if args.max_cells and adata.n_obs > args.max_cells:
            console.print(f"  Subsampling to {args.max_cells} cells")
            adata = adata[:args.max_cells, :]
            cell_metadata = adata.obs

        # Create DataLoader
        console.print("Creating DataLoader...")
        dataloader = tahoe_adapter.create_dataloader_from_adata(
            adata,
            batch_size=EXTRACTION_BATCH_SIZE,
            max_length=EXTRACTION_MAX_LENGTH,
            gene_id_key=args.gene_id_key,
            num_workers=EXTRACTION_NUM_WORKERS,
        )

        # Extract activations
        console.print("Extracting activations from Tahoe X1...")
        activations, gene_ids = tahoe_adapter.extract_activations(dataloader)
        console.print(f"✓ Extracted activations: {activations.shape}")

        # Save for future use
        torch.save(
            {
                "activations": activations,
                "gene_ids": gene_ids,
            },
            activations_path,
        )
        console.print(f"✓ Saved activations to {activations_path}")

    # Step 4: Extract gene-level associations
    console.print("\n[bold]Step 4: Extracting gene-level feature associations[/bold]")
    console.print("Analyzing which genes activate each SAE feature...")

    gene_results = analyzer.extract_feature_gene_associations(
        sae=trainer.sae,
        activations=activations,
        gene_ids=gene_ids,
        expressions=None,  # Could load from adata if available
        top_k=args.top_k_genes,
        device=args.device,
    )

    console.print(f"✓ Analyzed {trainer.d_sae} features")

    # Print some examples
    console.print("\n[bold]Sample Results:[/bold]")
    for feat_idx in [0, 1, 2]:
        if feat_idx in gene_results and gene_results[feat_idx]["top_genes"]:
            top_genes = gene_results[feat_idx]["top_genes"][:5]
            n_act = gene_results[feat_idx]["n_activations"]
            console.print(f"\nFeature {feat_idx} (activated {n_act} times):")
            console.print(f"  Top genes:")
            for gene, score in top_genes:
                console.print(f"    - {gene}: {score:.4f}")

    # Save gene results
    gene_results_path = args.output / "feature_gene_associations.json"
    analyzer.save_results(gene_results, gene_results_path)
    console.print(f"\n✓ Saved gene associations to {gene_results_path}")

    # Step 5: Extract cell-level associations
    console.print("\n[bold]Step 5: Extracting cell-level feature associations[/bold]")
    console.print("Analyzing which cells activate each SAE feature...")

    # Use cell metadata loaded in step 3
    if cell_metadata is not None:
        console.print(f"  Using cell metadata with {len(cell_metadata.columns)} columns")
    else:
        console.print("  [yellow]No cell metadata available[/yellow]")

    cell_results = analyzer.extract_top_activating_cells(
        sae=trainer.sae,
        activations=activations,
        cell_metadata=cell_metadata,
        top_k_cells=args.top_k_cells,
        device=args.device,
    )

    console.print(f"✓ Analyzed cell associations")

    # Print some examples
    if cell_metadata is not None and "cell_type" in cell_metadata.columns:
        console.print("\n[bold]Sample Cell Type Enrichment:[/bold]")
        for feat_idx in [0, 1, 2]:
            if feat_idx in cell_results and cell_results[feat_idx]["cell_types"]:
                cell_types = cell_results[feat_idx]["cell_types"]
                console.print(f"\nFeature {feat_idx} enriched in:")
                for cell_type, count in sorted(
                    cell_types.items(), key=lambda x: x[1], reverse=True
                )[:3]:
                    console.print(f"  - {cell_type}: {count} cells")

    # Save cell results
    cell_results_path = args.output / "feature_cell_associations.json"
    analyzer.save_results(cell_results, cell_results_path)
    console.print(f"\n✓ Saved cell associations to {cell_results_path}")

    # Step 6: Create summary
    console.print("\n[bold]Step 6: Creating summary[/bold]")

    summary = {
        "sae_checkpoint": str(args.sae),
        "data_file": str(args.data),
        "n_features": trainer.d_sae,
        "n_cells_analyzed": activations.shape[0],
        "top_k_genes": args.top_k_genes,
        "top_k_cells": args.top_k_cells,
        "active_features": sum(
            1
            for results in gene_results.values()
            if results["n_activations"] > 0
        ),
    }

    summary_path = args.output / "extraction_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    console.print(f"✓ Saved summary to {summary_path}")

    # Final summary
    console.print("\n[bold green]Feature Extraction Complete![/bold green]")
    console.print(f"Results saved to: {args.output}")
    console.print(f"\nSummary:")
    console.print(f"  - Features analyzed: {summary['n_features']}")
    console.print(f"  - Active features: {summary['active_features']}")
    console.print(f"  - Cells analyzed: {summary['n_cells_analyzed']}")
    console.print(f"\nNext step:")
    console.print(f"  Create visualization dashboard:")
    console.print(f"     python scripts/03_create_dashboard.py --results {args.output}")


if __name__ == "__main__":
    main()
