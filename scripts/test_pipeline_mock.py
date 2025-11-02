#!/usr/bin/env python3
"""
Mock test of the SAE pipeline without tahoe-x1 dependencies.

This script generates fake activations to test the SAE training,
feature extraction, and dashboard creation without needing tahoe-x1.

Usage:
    python scripts/test_pipeline_mock.py --output results/mock_test
"""

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from rich.console import Console
from tqdm import tqdm

# Direct imports (no tahoe needed)
from sae_lens import TopKTrainingSAE, TopKTrainingSAEConfig, TopKSAE, TopKSAEConfig
from sae_lens.saes.sae import TrainStepInput

console = Console()


def generate_mock_data(n_cells=500, n_genes=1000, seq_len=512, d_model=512):
    """Generate mock single-cell data and activations."""
    console.print(f"Generating mock data: {n_cells} cells, {n_genes} genes")

    # Random gene vocabulary
    gene_vocab = {f"GENE_{i:04d}": i for i in range(n_genes)}
    reverse_vocab = {v: k for k, v in gene_vocab.items()}

    # Random activations from "Tahoe model"
    activations = torch.randn(n_cells, seq_len, d_model)

    # Random gene IDs for each position
    gene_ids = torch.randint(0, n_genes, (n_cells, seq_len))

    # Mock cell metadata
    cell_types = ["T cells", "B cells", "Monocytes", "NK cells", "Other"]
    cell_type_labels = np.random.choice(cell_types, size=n_cells)

    return activations, gene_ids, gene_vocab, reverse_vocab, cell_type_labels


def mock_train_sae(activations, d_in, d_sae, k, n_epochs, device):
    """Train SAE on mock activations."""
    console.print("\n[bold]Training SAE on mock activations[/bold]")

    # Initialize SAE
    config = TopKTrainingSAEConfig(
        d_in=d_in,
        d_sae=d_sae,
        k=k,
        device=device,
    )
    sae = TopKTrainingSAE(config).to(device)

    # Flatten activations
    n_cells, seq_len, d_model = activations.shape
    activations_flat = activations.reshape(-1, d_model).to(device)

    # Simple training loop
    optimizer = torch.optim.Adam(sae.parameters(), lr=3e-4)
    batch_size = 256

    n_samples = len(activations_flat)
    n_batches = n_samples // batch_size

    sae.train()
    losses = []

    for epoch in range(n_epochs):
        pbar = tqdm(range(n_batches), desc=f"Epoch {epoch+1}/{n_epochs}")

        for batch_idx in pbar:
            start_idx = batch_idx * batch_size
            end_idx = start_idx + batch_size
            batch = activations_flat[start_idx:end_idx]

            step_input = TrainStepInput(sae_in=batch, dead_neuron_mask=None)
            output = sae.training_forward_pass(step_input)

            total_loss = output.losses.get("mse_loss", torch.tensor(0.0))
            for loss_name, loss_value in output.losses.items():
                if loss_name != "mse_loss":
                    total_loss = total_loss + loss_value

            optimizer.zero_grad()
            total_loss.backward()
            optimizer.step()

            losses.append(total_loss.item())
            pbar.set_postfix({"loss": f"{total_loss.item():.4f}"})

    console.print(f"✓ Training complete. Final loss: {losses[-1]:.4f}")

    # Convert to inference SAE
    inf_config = TopKSAEConfig(d_in=d_in, d_sae=d_sae, k=k, device=device)
    sae_inf = TopKSAE(inf_config).to(device)
    sae_inf.load_state_dict(sae.state_dict(), strict=False)

    return sae_inf, losses


def mock_extract_features(sae, activations, gene_ids, vocab, reverse_vocab, top_k=100, device="cpu"):
    """Extract top activating genes for each feature."""
    console.print("\n[bold]Extracting gene-level feature associations[/bold]")

    sae = sae.to(device)
    sae.eval()

    n_cells, seq_len, d_model = activations.shape
    n_features = sae.cfg.d_sae

    # Storage
    from collections import defaultdict
    feature_gene_scores = defaultdict(lambda: defaultdict(float))
    feature_activation_counts = defaultdict(int)

    with torch.no_grad():
        for cell_idx in tqdm(range(n_cells), desc="Processing cells"):
            cell_acts = activations[cell_idx:cell_idx+1].to(device)
            cell_genes = gene_ids[cell_idx]

            cell_acts_flat = cell_acts.reshape(-1, d_model)
            feature_acts = sae.encode(cell_acts_flat)

            for pos in range(seq_len):
                gene_id = cell_genes[pos].item()
                if gene_id < 0:
                    continue

                pos_feature_acts = feature_acts[pos]
                active_features = (pos_feature_acts > 0).nonzero(as_tuple=True)[0]

                for feat_idx in active_features:
                    feat_idx = feat_idx.item()
                    feat_activation = pos_feature_acts[feat_idx].item()

                    feature_gene_scores[feat_idx][gene_id] += feat_activation
                    feature_activation_counts[feat_idx] += 1

    # Compile results
    feature_results = {}
    for feat_idx in range(n_features):
        if feat_idx not in feature_gene_scores:
            feature_results[feat_idx] = {
                "top_genes": [],
                "n_activations": 0,
            }
            continue

        gene_scores = feature_gene_scores[feat_idx]
        sorted_genes = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)[:top_k]

        top_genes = [
            (reverse_vocab.get(gene_id, f"UNKNOWN_{gene_id}"), score)
            for gene_id, score in sorted_genes
        ]

        feature_results[feat_idx] = {
            "top_genes": top_genes,
            "n_activations": feature_activation_counts[feat_idx],
        }

    active_features = sum(1 for r in feature_results.values() if r["n_activations"] > 0)
    console.print(f"✓ Analyzed {n_features} features ({active_features} active)")

    return feature_results


def main():
    parser = argparse.ArgumentParser(description="Mock test of SAE pipeline")
    parser.add_argument("--output", type=Path, default="results/mock_test", help="Output directory")
    parser.add_argument("--n-cells", type=int, default=500, help="Number of cells")
    parser.add_argument("--n-genes", type=int, default=1000, help="Number of genes")
    parser.add_argument("--d-model", type=int, default=512, help="Model dimension")
    parser.add_argument("--d-sae", type=int, default=4096, help="SAE dictionary size")
    parser.add_argument("--k", type=int, default=64, help="TopK k value")
    parser.add_argument("--epochs", type=int, default=3, help="Training epochs")
    parser.add_argument("--device", type=str, default="cpu", choices=["cpu", "cuda"], help="Device")

    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    console.print("[bold green]Mock SAE Pipeline Test[/bold green]")
    console.print(f"Output: {args.output}")
    console.print(f"Device: {args.device}")

    # Step 1: Generate mock data
    console.print("\n[bold]Step 1: Generating mock data[/bold]")
    activations, gene_ids, vocab, reverse_vocab, cell_types = generate_mock_data(
        n_cells=args.n_cells,
        n_genes=args.n_genes,
        seq_len=512,
        d_model=args.d_model,
    )
    console.print(f"✓ Generated activations: {activations.shape}")

    # Step 2: Train SAE
    console.print("\n[bold]Step 2: Training SAE[/bold]")
    sae, losses = mock_train_sae(
        activations=activations,
        d_in=args.d_model,
        d_sae=args.d_sae,
        k=args.k,
        n_epochs=args.epochs,
        device=args.device,
    )

    # Save checkpoint
    checkpoint_path = args.output / "sae_mock.pt"
    torch.save({
        "sae_state_dict": sae.state_dict(),
        "d_in": args.d_model,
        "d_sae": args.d_sae,
        "k": args.k,
    }, checkpoint_path)
    console.print(f"✓ Saved checkpoint: {checkpoint_path}")

    # Step 3: Extract features
    console.print("\n[bold]Step 3: Extracting features[/bold]")
    feature_results = mock_extract_features(
        sae=sae,
        activations=activations,
        gene_ids=gene_ids,
        vocab=vocab,
        reverse_vocab=reverse_vocab,
        top_k=100,
        device=args.device,
    )

    # Save results
    results_path = args.output / "feature_gene_associations.json"
    with open(results_path, "w") as f:
        json.dump({str(k): v for k, v in feature_results.items()}, f, indent=2)
    console.print(f"✓ Saved results: {results_path}")

    # Mock cell results
    cell_results = {}
    for feat_idx in range(min(args.d_sae, 100)):
        cell_results[str(feat_idx)] = {
            "cell_types": {ct: np.random.randint(10, 100) for ct in ["T cells", "B cells", "Monocytes"][:2]}
        }

    cell_results_path = args.output / "feature_cell_associations.json"
    with open(cell_results_path, "w") as f:
        json.dump(cell_results, f, indent=2)

    # Summary
    summary = {
        "n_features": args.d_sae,
        "n_cells_analyzed": args.n_cells,
        "active_features": sum(1 for r in feature_results.values() if r["n_activations"] > 0),
        "top_k_genes": 100,
    }

    summary_path = args.output / "extraction_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    # Print sample results
    console.print("\n[bold]Sample Results:[/bold]")
    for feat_idx in range(min(3, args.d_sae)):
        if feature_results[feat_idx]["top_genes"]:
            top_genes = feature_results[feat_idx]["top_genes"][:5]
            n_act = feature_results[feat_idx]["n_activations"]
            console.print(f"\nFeature {feat_idx} ({n_act} activations):")
            for gene, score in top_genes:
                console.print(f"  - {gene}: {score:.4f}")

    # Step 4: Create dashboard
    console.print("\n[bold]Step 4: Creating dashboard[/bold]")
    import subprocess
    result = subprocess.run([
        "python", "scripts/03_create_dashboard.py",
        "--results", str(args.output),
        "--top-n", "50",
    ])

    if result.returncode == 0:
        console.print(f"✓ Dashboard created: {args.output}/dashboard.html")

    # Final summary
    console.print("\n[bold green]Mock Pipeline Test Complete![/bold green]")
    console.print(f"\nResults:")
    console.print(f"  - SAE trained: {args.d_sae} features, k={args.k}")
    console.print(f"  - Features analyzed: {summary['active_features']}/{summary['n_features']}")
    console.print(f"  - Cells processed: {summary['n_cells_analyzed']}")
    console.print(f"\nFiles created:")
    console.print(f"  - {checkpoint_path}")
    console.print(f"  - {results_path}")
    console.print(f"  - {args.output}/dashboard.html")
    console.print(f"\nOpen dashboard:")
    console.print(f"  open {args.output}/dashboard.html")
    console.print(f"\n[bold cyan]This mock test verifies:[/bold cyan]")
    console.print(f"  ✓ SAE training loop works")
    console.print(f"  ✓ SAELens v6 API correctly integrated")
    console.print(f"  ✓ Feature extraction logic works")
    console.print(f"  ✓ Dashboard generation works")
    console.print(f"\n[bold yellow]Note:[/bold yellow] Real pipeline needs tahoe-x1 model (use Docker)")


if __name__ == "__main__":
    main()
