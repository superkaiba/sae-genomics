#!/usr/bin/env python3
"""Quick script to check d_model for 1.3B Tahoe model."""

import sys
from pathlib import Path

# Add tahoe-x1 to path
sys.path.insert(0, str(Path('external/tahoe-x1')))

try:
    from tahoe_x1.model import ComposerTX

    print("Loading 1.3B model config from HuggingFace...")

    # Load model config (will download weights but we only need config)
    model, vocab, _, _ = ComposerTX.from_hf("tahoebio/Tahoe-x1", "1.3b")

    # Get config info
    d_model = model.model.d_model
    n_layers = model.model.n_layers
    n_heads = model.model.n_heads
    expansion_ratio = model.model.expansion_ratio

    print(f"\nTahoe X1 1.3B Model Architecture:")
    print(f"  d_model (residual stream): {d_model}")
    print(f"  n_layers: {n_layers}")
    print(f"  n_heads: {n_heads}")
    print(f"  expansion_ratio: {expansion_ratio}")
    print(f"  MLP dimension: {d_model * expansion_ratio}")

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
