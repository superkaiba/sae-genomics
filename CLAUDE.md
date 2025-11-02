# CLAUDE.md - AI Assistant Context

This file provides context for AI assistants working on this repository.

## Project Overview

**SAE Genomics** - Train sparse autoencoders (SAEs) on the tahoe x1 single-cell foundation model to extract interpretable biological features.

**Purpose**: Mechanistic interpretability for single-cell genomics models. SAEs decompose tahoe x1's learned representations into sparse, interpretable features that can be validated against biological databases.

## Repository Architecture

### High-Level Pipeline Flow

```
Single-cell data (.h5ad)
    ↓
Tahoe X1 Model (transformer, processes gene expression)
    ↓
Intermediate Layer Activations (MLP outputs at specific layers)
    ↓
Sparse Autoencoder (TopK with k=64, learns d_sae=16384 features)
    ↓
Feature Analysis (both gene-level and cell-level)
    ↓
Biological Validation (10 databases) + Interactive Dashboard
```

### Key Components

**1. Tahoe X1 Integration** (`src/sae_genomics/models/`)
- `tahoe_adapter.py`: Loads tahoe-x1 from HuggingFace, extracts activations
- `model_config.py`: Configuration for different model sizes (70M, 1.3B, 3B)
- Hook points: Extract at specific transformer layers (e.g., `blocks.0.hook_mlp_out`)
- Lazy imports to avoid circular dependencies with tahoe-x1

**2. SAE Training** (`src/sae_genomics/training/`)
- `trainer.py`: Uses SAELens v6 API (TopKTrainingSAE/TopKSAEConfig)
- `config.py`: Training hyperparameters and settings
- Training loop uses `training_forward_pass()` with `TrainStepInput`/`TrainStepOutput`
- Saves both training SAE and inference SAE checkpoints

**3. Feature Analysis** (`src/sae_genomics/utils/`)
- `feature_analysis.py`: Extracts gene and cell associations
- **Gene-level**: `score = feature_activation × gene_expression` (which genes activate each feature)
- **Cell-level**: Which cells/cell types most strongly activate each feature
- Results saved as JSON for downstream processing

**4. Visualization** (`scripts/03_create_dashboard.py`)
- Creates self-contained HTML dashboard
- Shows top genes per feature, cell type enrichment
- Search functionality, modern dark theme
- No external dependencies (all CSS/JS inline)

**5. Validation Framework** (`src/sae_genomics/validation/`)
- Planned: 10 biological databases for feature validation
- Tier 1: Open Targets, HPO, GWAS Catalog
- Tier 2: DisGeNET, Monarch, STRING
- Tier 3: OMIM, Orphanet, ClinVar, GenCC
- **Status**: Structure in place, implementation pending

## Key Design Decisions

### Why TopK SAE?

- **Sparse activation**: Only top-k=64 features active per input
- **Interpretability**: Each feature should capture a distinct biological pattern
- **SAELens v6**: Uses TopKTrainingSAE for training, TopKSAE for inference
- **No L1 penalty needed**: TopK enforces sparsity directly

### Why Extract Both Gene & Cell Associations?

**Gene-level**: "Feature 42 activates for [IL6, TNF, CXCL8]" → biological process
**Cell-level**: "Feature 88 activates in T cells, NK cells" → cell type identity

Both perspectives are valuable:
- Genes → biological pathways, disease associations
- Cells → cell type markers, state identification

### Data Format

**Input**: AnnData files (`.h5ad`)
- `adata.X`: Gene expression matrix (cells × genes)
- `adata.var['ensembl_id']`: Gene identifiers
- `adata.obs['cell_type']`: Optional cell annotations

**Why AnnData?**: Standard in single-cell biology, used by scanpy, tahoe-x1

### Activation Extraction Strategy

Extract from **intermediate layers**, not final embeddings:
- Hook points: `blocks.{layer_num}.hook_mlp_out`
- Earlier layers: Low-level gene patterns
- Later layers: High-level cell state
- Default: Extract from multiple layers (0, 3, 6, 9, 11 for 70M model)

## File Structure Guide

### Scripts (User-Facing)

```
scripts/
├── download_test_data.py       # Gets PBMC 3k dataset
├── 01_train_sae.py            # Step 1: Train SAE
├── 02_extract_features.py     # Step 2: Extract gene/cell associations
├── 03_create_dashboard.py     # Step 3: Create visualization
├── run_pipeline.py            # Master: Runs all 3 steps
└── test_pipeline_mock.py      # Mock test without tahoe (for debugging)
```

### Source Code (Implementation)

```
src/sae_genomics/
├── models/
│   ├── tahoe_adapter.py       # Model loading, activation extraction
│   └── model_config.py        # Configuration dataclasses
├── training/
│   ├── trainer.py             # SAE training logic
│   └── config.py              # Training configuration
├── validation/                # Biological validation (planned)
│   ├── databases/            # 10 database clients (to implement)
│   ├── enrichment/           # Statistical tests (to implement)
│   └── visualization/        # Plots (to implement)
├── data/                     # Data loading (minimal, uses tahoe utilities)
└── utils/
    └── feature_analysis.py   # Gene/cell association extraction
```

### Configuration

```
configs/
├── models/
│   ├── tx1_70m.yaml          # 70M model (prototyping)
│   └── tx1_3b.yaml           # 3B model (production)
├── training/
│   └── default.yaml          # SAE training config
└── validation/
    └── databases.yaml        # Database API endpoints and settings
```

## Important Implementation Details

### SAELens v6 API

**Training:**
```python
config = TopKTrainingSAEConfig(d_in=512, d_sae=16384, k=64)
sae = TopKTrainingSAE(config)
step_input = TrainStepInput(sae_in=batch, dead_neuron_mask=None)
output = sae.training_forward_pass(step_input)
loss = output.losses["mse_loss"]
```

**Inference:**
```python
config = TopKSAEConfig(d_in=512, d_sae=16384, k=64)
sae = TopKSAE(config)
feature_acts = sae.encode(activations)  # (batch, d_sae)
reconstructed = sae.decode(feature_acts)  # (batch, d_in)
```

### Lazy Imports Pattern

Used throughout to avoid circular dependencies and early loading of heavy libraries:

```python
# In __init__.py
def __getattr__(name):
    if name == "TahoeModelAdapter":
        from sae_genomics.models.tahoe_adapter import TahoeModelAdapter
        return TahoeModelAdapter
    raise AttributeError(...)

# In tahoe_adapter.py
def _import_tahoe():
    """Lazy import of tahoe_x1 to avoid dependency issues."""
    from tahoe_x1.model import ComposerTX
    from tahoe_x1.utils.util import loader_from_adata
    return ComposerTX, loader_from_adata
```

### Git Submodules

```
external/
├── tahoe-x1/      # Tahoe model implementation
└── SAELens/       # SAE training framework
```

**Important**: Use `git clone --recursive` or run `git submodule update --init --recursive`

## Common Tasks

### Adding a New Script

1. Create in `scripts/`
2. Make executable: `chmod +x scripts/new_script.py`
3. Add shebang: `#!/usr/bin/env python3`
4. Import from `sae_genomics` package
5. Use `rich.console` for output
6. Add argparse for CLI

### Adding a New Module

1. Create in appropriate `src/sae_genomics/` subdirectory
2. Add to `__init__.py` (use lazy imports if heavy dependencies)
3. Add type hints
4. Write docstrings
5. Add to `__all__`

### Implementing Database Clients

Template for `src/sae_genomics/validation/databases/new_db.py`:

```python
from typing import Dict, List
import requests

class NewDBClient:
    def __init__(self, api_url: str, cache_dir: str):
        self.api_url = api_url
        self.cache_dir = cache_dir

    def query_genes(self, gene_list: List[str]) -> Dict:
        # Query API
        # Cache results
        # Return associations
        pass
```

## Testing Strategy

### Local Testing (macOS)

**Blocked**: tahoe-x1 requires CUDA/Linux
- Complex dependency chain: llm-foundry → flash-attn → CUDA
- Use mock tests or skip to RunPod

### RunPod/GPU Testing

**Repository**: https://github.com/superkaiba/sae-genomics

**Clone on RunPod:**
```bash
git clone --recursive https://github.com/superkaiba/sae-genomics.git
```

**Quick test:**
```bash
cd sae-genomics
uv venv && source .venv/bin/activate
uv pip install scanpy && uv pip install -e .
python scripts/download_test_data.py
python scripts/run_pipeline.py \
  --data data/pbmc3k_test.h5ad \
  --output results/test \
  --steps 500 \
  --max-cells 500 \
  --device cuda
```

### Mock Testing (No Tahoe)

Use `scripts/test_pipeline_mock.py` to test SAE training and feature extraction logic without tahoe dependencies.

## Dependency Management

### Core Dependencies

```
torch>=2.1.0           # Deep learning
transformers>=4.35.0   # HuggingFace models
sae-lens>=6.20.0       # SAE training
scanpy>=1.9.0          # Single-cell analysis
```

### Tahoe-x1 Dependencies

```
llm-foundry[gpu]>=0.17.1   # LLM training framework
composer>=0.32.0           # Training infrastructure
mosaicml-streaming>=0.7.0  # Data streaming
omegaconf>=2.3.0          # Configuration
boto3>=1.40.0             # AWS S3 access
```

### Installation Order (Important!)

1. `uv pip install scanpy` (separate to avoid conflicts)
2. `uv pip install -e .` (installs rest from pyproject.toml)

**Note**: On macOS, llm-foundry installation may fail due to sentencepiece/flash-attn. This is expected - use Docker or Linux.

## Known Issues & Workarounds

### 1. macOS Segfaults

**Issue**: tahoe model loading causes exit code 139
**Cause**: flash-attn compiled for CUDA, not compatible with macOS
**Solution**: Use Docker or Linux+GPU

### 2. Sentencepiece Build Failures

**Issue**: llm-foundry requires sentencepiece which needs cmake
**Solution**: `brew install cmake` OR use prebuilt wheels: `pip install sentencepiece`

### 3. Import Errors on Package Load

**Issue**: Circular imports or missing modules on `import sae_genomics`
**Solution**: All heavy imports are lazy-loaded via `__getattr__`

### 4. SAELens API Changes

**Version**: We use SAELens v6.20.1
**Key changes from older versions**:
- Use `TopKTrainingSAE` not `SAE`
- Use `training_forward_pass()` not `forward()`
- Use `.encode()` not `(sae_out, features, loss)`
- Config classes: `TopKSAEConfig`, `TopKTrainingSAEConfig`

## Performance Notes

### Model Sizes

| Model | Params | Context | d_model | GPU Memory | Speed |
|-------|--------|---------|---------|------------|-------|
| 70M | 70M | 1024 | 512 | ~4GB | Fast (prototyping) |
| 1.3B | 1.3B | 2048 | 2048 | ~8GB | Medium |
| 3B | 3B | 2048 | 2048 | ~16GB | Slow (production) |

### SAE Training

- SAE adds ~2-4GB GPU memory
- TopK is memory-efficient (sparse activations)
- Typical training: 10K-50K steps
- Time: ~30-60 min for 10K steps on H100

### Data Processing

- Extract activations in batches
- Cache to disk to avoid recomputation
- Saved as `results/{experiment}/activations.pt`
- Reusable across feature extraction runs

## Code Conventions

### Imports

- Standard library first
- Third-party packages second
- Local imports last
- Lazy imports for heavy dependencies

### Type Hints

- Use `from __future__ import annotations`
- Use `TYPE_CHECKING` for circular import type hints
- All functions should have type hints
- Use jaxtyping for tensor shapes when relevant

### Configuration

- YAML for all configuration files
- Load with `yaml.safe_load()`
- Use dataclasses for config objects
- Allow CLI overrides

### Error Handling

- Use `rich.console` for user-facing output
- Log errors with context
- Fail fast with clear error messages
- Suggest solutions in error messages

## Future Work (Planned)

### Immediate Next Steps

1. **Implement database clients** (`src/sae_genomics/validation/databases/`)
   - Open Targets, HPO, GWAS Catalog (Tier 1)
   - DisGeNET, Monarch, STRING (Tier 2)
   - OMIM, Orphanet, ClinVar, GenCC (Tier 3)

2. **Enrichment analysis** (`src/sae_genomics/validation/enrichment/`)
   - Fisher's exact test for disease/pathway enrichment
   - Multiple testing correction (FDR)
   - HPO hierarchy-aware phenotype enrichment

3. **Biological validation pipeline**
   - Integrate database queries
   - Generate validation reports
   - Visualize enrichment results

### Advanced Features

1. **Distributed training** for 3B model
2. **Feature steering**: Modify SAE features, observe effect
3. **Cross-model comparison**: Compare features across model sizes
4. **Perturbation analysis**: Use tahoe's perturbation data
5. **Web dashboard**: Replace static HTML with interactive app (Streamlit)

## Quick Reference

### File Locations

**Main scripts**: `scripts/01_train_sae.py`, `02_extract_features.py`, `03_create_dashboard.py`
**Orchestration**: `scripts/run_pipeline.py`
**Core logic**: `src/sae_genomics/{models,training,utils}/`
**Configs**: `configs/{models,training,validation}/`
**Test data**: `data/pbmc3k_test.h5ad` (download with `scripts/download_test_data.py`)

### Key Commands

**Full pipeline:**
```bash
python scripts/run_pipeline.py --data data.h5ad --output results/exp1 --steps 10000
```

**Individual steps:**
```bash
python scripts/01_train_sae.py --data data.h5ad --output results/exp1 --steps 10000
python scripts/02_extract_features.py --sae results/exp1/sae_final.pt --data data.h5ad
python scripts/03_create_dashboard.py --results results/exp1
```

**Quick test:**
```bash
python scripts/run_pipeline.py --data data.h5ad --output results/test --steps 500 --max-cells 500
```

### Configuration Override

CLI arguments override config files:
```bash
--config configs/training/custom.yaml  # Use custom config
--model-size 70m                       # Override model
--steps 5000                           # Override training steps
--device cuda                          # Override device
```

## Data Flow Details

### Tahoe X1 Processing

1. **Input**: AnnData with gene expression (cells × genes)
2. **Tokenization**: Genes → vocabulary IDs
3. **Model forward**: Through transformer layers
4. **Hook extraction**: Capture MLP outputs at specified layers
5. **Output**: Activations tensor (cells, seq_len, d_model)

### SAE Training

1. **Input**: Activations (cells × seq_len × d_model)
2. **Flatten**: Reshape to (cells * seq_len, d_model)
3. **Training**: TopK SAE learns sparse decomposition
4. **Output**: Checkpoint with W_enc, W_dec, b_enc, b_dec weights

### Feature Extraction

1. **Input**: Trained SAE + activations + gene IDs
2. **For each cell**:
   - Run through SAE.encode() → feature activations
   - For each position (gene) in sequence:
     - Track which features activate
     - Weight by feature_activation × gene_expression
3. **Aggregate**: Sum scores across all cells per (feature, gene) pair
4. **Output**: Top-k genes for each feature

### Dashboard Creation

1. **Input**: Feature-gene associations + feature-cell associations (JSON)
2. **Processing**: Sort features by activation frequency
3. **HTML generation**: Create self-contained HTML
4. **Output**: Interactive dashboard with search

## Environment Requirements

### Development (macOS)

- Python 3.10+
- uv package manager
- Can edit code, write tests
- **Cannot run full pipeline** (needs GPU/Linux)

### Production (Linux + GPU)

- NVIDIA GPU (A100, H100, or similar)
- CUDA 12.1+
- Python 3.10-3.11
- Docker (recommended) or native install

### Testing Environments

- **RunPod**: GPU cloud instances
- **Google Colab**: Free GPU tier
- **AWS/GCP/Azure**: Cloud GPU instances
- **Local workstation**: If you have NVIDIA GPU

## Debugging Guide

### Common Errors

**"ModuleNotFoundError: No module named 'tahoe_x1'"**
- Submodules not initialized
- Run: `git submodule update --init --recursive`

**"ModuleNotFoundError: No module named 'llmfoundry'"**
- llm-foundry not installed
- Run: `uv pip install llm-foundry` (may need cmake)

**Exit code 139 (Segfault)**
- Native library crash, usually on macOS
- Solution: Use Docker or Linux+GPU

**"CUDA out of memory"**
- Reduce `--max-cells` or `--batch-size`
- Use smaller model (70M instead of 3B)
- Use gradient accumulation

**"No genes matched to vocabulary"**
- Gene IDs don't match tahoe vocabulary
- Check `--gene-id-key` argument
- Tahoe uses Ensembl IDs by default

### Logging Locations

- Console output: Rich formatted progress
- Training metrics: `results/{exp}/training_metrics.json`
- Extraction summary: `results/{exp}/extraction_summary.json`
- Model info: `results/{exp}/{model,train}_config.json`

## Performance Optimization

### For Fast Iteration

- Use 70M model
- Limit cells: `--max-cells 1000`
- Short training: `--steps 1000`
- CPU acceptable for small tests

### For Production

- Use 3B model
- Full dataset
- Long training: `--steps 20000-50000`
- Multiple GPUs if available
- Save activation cache for reuse

### Memory Management

- Activations cached to disk after extraction
- Reuse with `--use-saved-activations` flag
- Clear results/ directory if running low on space
- Use `--max-cells` to limit memory usage

## Contributing Guidelines

### Code Style

- Black formatting (line length 100)
- Ruff linting
- Type hints required
- Docstrings for all public functions

### Commit Messages

- Descriptive, present tense
- Include context and reasoning
- Reference issues if applicable
- Use conventional commits format

### Testing

- Add tests in `tests/` directory
- Use pytest
- Test both with and without tahoe (mock tests)
- CI/CD via GitHub Actions (future)

## Resources & References

### Documentation

- **This project**: See README.md, docs/
- **Tahoe X1**: https://github.com/tahoebio/tahoe-x1
- **SAELens**: https://github.com/jbloomAus/SAELens
- **Paper**: Tahoe-x1 bioRxiv preprint

### External Links

- **GitHub**: https://github.com/superkaiba/sae-genomics
- **HuggingFace**: https://huggingface.co/tahoebio/Tahoe-x1 (model weights)
- **CellxGene**: https://cellxgene.cziscience.com/ (data source)

### Key Papers

1. "Sparse Autoencoders Find Interpretable Features in Language Models" (Anthropic)
2. "Scaling Monosemanticity" (Anthropic)
3. "Sparse Autoencoders Reveal Interpretable Features in Single-Cell Foundation Models" (bioRxiv 2024)
4. "Tahoe-x1: Scaling Perturbation-Trained Single-Cell Foundation Models" (bioRxiv 2025)

## Version Information

- **Project version**: 0.1.0
- **SAELens**: v6.20.1
- **Tahoe X1**: v0.1.0-17
- **Python**: 3.10-3.13 supported
- **Last updated**: 2025-11-02

## Contact & Support

- **GitHub Issues**: https://github.com/superkaiba/sae-genomics/issues
- **Tahoe X1 Issues**: https://github.com/tahoebio/tahoe-x1/issues
- **SAELens Issues**: https://github.com/jbloomAus/SAELens/issues

---

**For AI Assistants**: This repository is production-ready for the core pipeline (training, extraction, visualization). The validation framework (database clients, enrichment analysis) is structurally planned but not yet implemented. Priority should be given to testing the existing pipeline on GPU hardware, then implementing the biological validation components.
