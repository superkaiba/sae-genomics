"""Feature analysis utilities for SAE interpretation."""

from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm


class FeatureAnalyzer:
    """Analyze SAE features to find top activating genes."""

    def __init__(self, vocab: Dict[str, int], reverse_vocab: Dict[int, str] = None):
        """Initialize analyzer.

        Args:
            vocab: Gene name to ID mapping
            reverse_vocab: ID to gene name mapping (auto-created if None)
        """
        self.vocab = vocab
        if reverse_vocab is None:
            self.reverse_vocab = {v: k for k, v in vocab.items()}
        else:
            self.reverse_vocab = reverse_vocab

    def extract_feature_gene_associations(
        self,
        sae,
        activations: torch.Tensor,
        gene_ids: torch.Tensor,
        expressions: torch.Tensor = None,
        top_k: int = 100,
        device: str = "cuda" if torch.cuda.is_available() else "cpu",
    ) -> Dict[int, Dict]:
        """Extract gene associations for each SAE feature.

        Strategy:
        1. For each cell, run through SAE to get feature activations
        2. For each feature, track which genes are present when it activates strongly
        3. Weight genes by both expression level and feature activation strength

        Args:
            sae: Trained SAE model
            activations: Tahoe activations (n_cells, seq_len, d_model)
            gene_ids: Gene IDs for each position (n_cells, seq_len)
            expressions: Gene expression values (n_cells, seq_len), optional
            top_k: Number of top genes to return per feature
            device: Device to run on

        Returns:
            Dictionary mapping feature_id -> {
                'top_genes': [(gene_name, score), ...],
                'gene_scores': {gene_name: score, ...},
                'n_activations': int,
            }
        """
        sae = sae.to(device)
        sae.eval()

        n_cells, seq_len, d_model = activations.shape
        n_features = sae.d_sae

        # Storage for gene-feature associations
        # feature_gene_scores[feature_id][gene_id] = cumulative score
        feature_gene_scores = defaultdict(lambda: defaultdict(float))
        feature_activation_counts = defaultdict(int)

        print(f"Analyzing {n_cells} cells across {n_features} features...")

        with torch.no_grad():
            for cell_idx in tqdm(range(n_cells), desc="Processing cells"):
                # Get activations for this cell
                cell_acts = activations[cell_idx: cell_idx + 1].to(device)  # (1, seq_len, d_model)
                cell_genes = gene_ids[cell_idx]  # (seq_len,)

                # Get expression values if provided
                if expressions is not None:
                    cell_expr = expressions[cell_idx]  # (seq_len,)
                else:
                    cell_expr = torch.ones_like(cell_genes, dtype=torch.float32)

                # Flatten for SAE
                cell_acts_flat = cell_acts.reshape(-1, d_model)  # (seq_len, d_model)

                # Run through SAE to get feature activations
                _, feature_acts, _ = sae(cell_acts_flat, return_loss=False)  # (seq_len, n_features)

                # For each position in the sequence
                for pos in range(seq_len):
                    gene_id = cell_genes[pos].item()
                    if gene_id < 0:  # Skip padding
                        continue

                    gene_expr = cell_expr[pos].item() if expressions is not None else 1.0
                    pos_feature_acts = feature_acts[pos]  # (n_features,)

                    # Find active features at this position
                    # (features with non-zero activation)
                    active_features = (pos_feature_acts > 0).nonzero(as_tuple=True)[0]

                    for feat_idx in active_features:
                        feat_idx = feat_idx.item()
                        feat_activation = pos_feature_acts[feat_idx].item()

                        # Score = feature_activation * gene_expression
                        # This captures: "how much does this gene contribute to this feature"
                        score = feat_activation * gene_expr

                        feature_gene_scores[feat_idx][gene_id] += score
                        feature_activation_counts[feat_idx] += 1

        # Convert to readable format
        feature_results = {}

        for feat_idx in tqdm(range(n_features), desc="Compiling results"):
            if feat_idx not in feature_gene_scores:
                feature_results[feat_idx] = {
                    "top_genes": [],
                    "gene_scores": {},
                    "n_activations": 0,
                }
                continue

            gene_scores = feature_gene_scores[feat_idx]

            # Sort genes by score
            sorted_genes = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)[
                :top_k
            ]

            # Convert to gene names
            top_genes = [
                (self.reverse_vocab.get(gene_id, f"UNKNOWN_{gene_id}"), score)
                for gene_id, score in sorted_genes
            ]

            # Also store all scores as gene names
            gene_scores_named = {
                self.reverse_vocab.get(gene_id, f"UNKNOWN_{gene_id}"): score
                for gene_id, score in gene_scores.items()
            }

            feature_results[feat_idx] = {
                "top_genes": top_genes,
                "gene_scores": gene_scores_named,
                "n_activations": feature_activation_counts[feat_idx],
            }

        return feature_results

    def extract_top_activating_cells(
        self,
        sae,
        activations: torch.Tensor,
        cell_metadata: pd.DataFrame = None,
        top_k_cells: int = 100,
        device: str = "cuda" if torch.cuda.is_available() else "cpu",
    ) -> Dict[int, Dict]:
        """Extract top activating cells for each feature.

        Args:
            sae: Trained SAE model
            activations: Tahoe activations (n_cells, seq_len, d_model)
            cell_metadata: DataFrame with cell annotations
            top_k_cells: Number of top cells to return per feature
            device: Device to run on

        Returns:
            Dictionary mapping feature_id -> {
                'top_cells': [(cell_idx, score), ...],
                'cell_types': Counter of cell types (if metadata provided),
            }
        """
        sae = sae.to(device)
        sae.eval()

        n_cells, seq_len, d_model = activations.shape
        n_features = sae.d_sae

        # Storage for cell-feature associations
        # feature_cell_scores[feature_id][cell_idx] = max activation
        feature_cell_scores = defaultdict(lambda: defaultdict(float))

        print(f"Analyzing top activating cells across {n_features} features...")

        with torch.no_grad():
            for cell_idx in tqdm(range(n_cells), desc="Processing cells"):
                cell_acts = activations[cell_idx: cell_idx + 1].to(device)
                cell_acts_flat = cell_acts.reshape(-1, d_model)

                # Get feature activations
                _, feature_acts, _ = sae(cell_acts_flat, return_loss=False)

                # Aggregate across sequence (max activation)
                max_feature_acts = feature_acts.max(dim=0).values  # (n_features,)

                # Store max activation for each feature
                for feat_idx in range(n_features):
                    activation_value = max_feature_acts[feat_idx].item()
                    if activation_value > 0:
                        feature_cell_scores[feat_idx][cell_idx] = activation_value

        # Compile results
        feature_results = {}

        for feat_idx in tqdm(range(n_features), desc="Compiling cell results"):
            if feat_idx not in feature_cell_scores:
                feature_results[feat_idx] = {
                    "top_cells": [],
                    "cell_types": {},
                }
                continue

            cell_scores = feature_cell_scores[feat_idx]

            # Sort cells by activation
            sorted_cells = sorted(cell_scores.items(), key=lambda x: x[1], reverse=True)[
                :top_k_cells
            ]

            # Get cell type distribution if metadata provided
            cell_types = {}
            if cell_metadata is not None and "cell_type" in cell_metadata.columns:
                from collections import Counter

                top_cell_indices = [c[0] for c in sorted_cells]
                top_cell_types = cell_metadata.iloc[top_cell_indices]["cell_type"].values
                cell_types = dict(Counter(top_cell_types))

            feature_results[feat_idx] = {
                "top_cells": sorted_cells,
                "cell_types": cell_types,
            }

        return feature_results

    @staticmethod
    def save_results(results: Dict, output_path: str):
        """Save feature analysis results."""
        import json

        # Convert to JSON-serializable format
        results_serializable = {}
        for feat_idx, feat_data in results.items():
            results_serializable[str(feat_idx)] = {
                "top_genes": feat_data.get("top_genes", []),
                "n_activations": feat_data.get("n_activations", 0),
                "top_cells": feat_data.get("top_cells", [])[:100],  # Limit to save space
                "cell_types": feat_data.get("cell_types", {}),
            }

        with open(output_path, "w") as f:
            json.dump(results_serializable, f, indent=2)

        print(f"Saved results to {output_path}")

    @staticmethod
    def load_results(input_path: str) -> Dict:
        """Load feature analysis results."""
        import json

        with open(input_path) as f:
            results = json.load(f)

        # Convert keys back to integers
        return {int(k): v for k, v in results.items()}
