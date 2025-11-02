# SAE Genomics

Training sparse autoencoders (SAEs) on the **tahoe x1** single-cell foundation model to extract interpretable biological features.

## Overview

This repository provides infrastructure for:
- Training sparse autoencoders on tahoe x1 model activations
- Comprehensive biological validation using 10 gene-disease/trait databases
- Feature interpretation and visualization
- Support for models from 70M to 3B parameters

## Architecture

### Tahoe X1 Model
[Tahoe X1](https://github.com/tahoebio/tahoe-x1) is a perturbation-trained single-cell foundation model:
- Processes transcriptomic data (266M+ cells)
- Transformer-based architecture (70M, 1.3B, 3B parameters)
- Generates cell and gene embeddings

### SAE Training
Using [SAELens](https://github.com/jbloomAus/SAELens), a comprehensive framework for training and analyzing sparse autoencoders.

### Biological Validation
Validates SAE features against 10 databases across 3 tiers:

**Tier 1 (Core):**
- Open Targets Platform - comprehensive disease associations
- Human Phenotype Ontology (HPO) - phenotype annotations
- GWAS Catalog - complex trait associations

**Tier 2 (Comprehensive):**
- DisGeNET - integrated disease-gene associations
- Monarch Initiative - cross-species integration
- STRING - protein interaction networks

**Tier 3 (Specialized):**
- OMIM - Mendelian diseases
- Orphanet - rare diseases
- ClinVar - clinical variant interpretations
- GenCC - gene-disease validity classifications

## Installation

### Prerequisites
- Python 3.10+
- [uv](https://github.com/astral-sh/uv) package manager
- CUDA-capable GPU (recommended for 3B model)

### Setup

1. **Clone the repository with submodules:**
```bash
git clone --recursive https://github.com/yourusername/sae-genomics.git
cd sae-genomics
```

If you already cloned without `--recursive`:
```bash
git submodule update --init --recursive
```

2. **Install dependencies with uv:**
```bash
# Install uv if you don't have it
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create virtual environment and install dependencies
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install the package in editable mode with all dependencies
uv pip install -e ".[all]"
```

3. **Set up environment variables:**
```bash
cp .env.example .env
# Edit .env with your API keys
```

Required API keys:
- `WANDB_API_KEY` - for experiment tracking
- `DISGENET_API_KEY` - for DisGeNET database access
- `OMIM_API_KEY` - for OMIM database access

## Quick Start

The complete pipeline consists of 3 steps:
1. Train SAE on Tahoe X1 activations
2. Extract top activating genes/cells for each feature
3. Create interactive visualization dashboard

### Option 1: Run Complete Pipeline (Recommended)

Run everything with one command:

```bash
python scripts/run_pipeline.py \
    --data path/to/your/data.h5ad \
    --model-size 70m \
    --output results/my_experiment \
    --steps 10000
```

**For quick testing:**
```bash
python scripts/run_pipeline.py \
    --data path/to/your/data.h5ad \
    --model-size 70m \
    --output results/test \
    --steps 1000 \
    --max-cells 1000
```

This will:
- ✓ Train SAE on Tahoe X1 activations
- ✓ Extract top genes and cells for each feature
- ✓ Create interactive HTML dashboard
- ✓ Save all results to `results/my_experiment/`

Open the dashboard: `open results/my_experiment/dashboard.html`

### Option 2: Run Steps Individually

**Step 1: Train SAE**
```bash
python scripts/01_train_sae.py \
    --data path/to/data.h5ad \
    --model-size 70m \
    --output results/my_experiment \
    --steps 10000
```

**Step 2: Extract Features**
```bash
python scripts/02_extract_features.py \
    --sae results/my_experiment/sae_final.pt \
    --data path/to/data.h5ad \
    --output results/my_experiment
```

**Step 3: Create Dashboard**
```bash
python scripts/03_create_dashboard.py \
    --results results/my_experiment
```

### What You Get

After running the pipeline, you'll have:

1. **Trained SAE** (`sae_final.pt`) - Sparse autoencoder with interpretable features
2. **Feature Analysis** (`feature_gene_associations.json`) - Top genes for each feature
3. **Cell Analysis** (`feature_cell_associations.json`) - Cell type enrichment per feature
4. **Interactive Dashboard** (`dashboard.html`) - Visualize all features

The dashboard shows:
- All SAE features ranked by activation frequency
- Top activating genes for each feature
- Cell type enrichment (if cell type annotations provided)
- Search functionality to filter by gene name

### Example with Your Data

Your data should be an AnnData file (`.h5ad`) with:
- `adata.X` - Gene expression matrix (cells × genes)
- `adata.var['ensembl_id']` - Ensembl gene IDs (or use `--gene-id-key` for different column)
- `adata.obs['cell_type']` - (Optional) Cell type annotations

```python
import scanpy as sc

# Load your data
adata = sc.read_h5ad("your_data.h5ad")
print(adata)  # Verify structure

# Run pipeline
!python scripts/run_pipeline.py \
    --data your_data.h5ad \
    --output results/exp1 \
    --steps 5000
```

## Pipeline Details

### How It Works

```
Single-cell data (h5ad) → Tahoe X1 → Layer Activations → SAE → Interpretable Features
                                            ↓
                                     (What we analyze)
```

**Step 1: Train SAE**
- Loads Tahoe X1 model from HuggingFace
- Processes your single-cell data through the model
- Extracts intermediate layer activations (not final embeddings)
- Trains a sparse autoencoder on these activations
- Saves trained SAE checkpoint

**Step 2: Extract Features**
- Loads trained SAE
- Runs cells through Tahoe X1 + SAE pipeline
- **Gene-level analysis**: Tracks which genes activate each SAE feature
  - For each feature, identifies top genes based on activation × expression
- **Cell-level analysis**: Identifies which cells strongly activate each feature
  - Analyzes cell type enrichment if annotations provided
- Saves feature-gene and feature-cell associations

**Step 3: Create Dashboard**
- Generates interactive HTML visualization
- Shows top genes for each feature
- Displays cell type enrichment
- Includes search functionality

### Understanding the Results

**SAE Features = Biological Patterns**

Each SAE feature captures a sparse pattern in the gene expression data. For example:

- **Feature 42** might activate for genes: IL6, TNF, CXCL8, CCL2
  - → Interpretation: "Inflammatory response feature"
- **Feature 88** might activate for genes: TOP2A, MKI67, PCNA, CDK1
  - → Interpretation: "Cell proliferation feature"
- **Feature 156** might activate in: T cells, NK cells
  - → Interpretation: "T cell identity feature"

**Top Activating Genes**

For each feature, the pipeline calculates:
```
score = feature_activation × gene_expression
```

This tells you: "which genes cause this feature to activate strongly"

### Command-Line Options

**run_pipeline.py options:**
```bash
--data PATH              # Path to h5ad file (required)
--output PATH            # Output directory (required)
--model-size {70m,1.3b,3b}  # Tahoe X1 size (default: 70m)
--steps INT              # Training steps (default: 10000)
--max-cells INT          # Limit cells for testing
--top-k-genes INT        # Top genes per feature (default: 100)
--top-k-cells INT        # Top cells per feature (default: 100)
--top-n-features INT     # Features in dashboard (default: 100)
--device {auto,cuda,cpu} # Compute device
--skip-training          # Skip training (use existing checkpoint)
--skip-extraction        # Skip feature extraction
--skip-dashboard         # Skip dashboard creation
```

**Example commands:**

```bash
# Quick test run (1000 cells, 1000 steps)
python scripts/run_pipeline.py \
    --data data.h5ad \
    --output results/test \
    --steps 1000 \
    --max-cells 1000

# Full run on 3B model
python scripts/run_pipeline.py \
    --data data.h5ad \
    --model-size 3b \
    --output results/large_model \
    --steps 50000

# Resume from existing checkpoint
python scripts/run_pipeline.py \
    --data data.h5ad \
    --output results/exp1 \
    --skip-training \
    --skip-extraction

# Extract more genes per feature
python scripts/02_extract_features.py \
    --sae results/exp1/sae_final.pt \
    --data data.h5ad \
    --top-k-genes 200
```

### Performance Tips

**For prototyping (fast iteration):**
- Use 70M model
- Limit to 1000-5000 cells (`--max-cells`)
- Train for 1000-5000 steps
- Takes ~10-30 minutes on GPU

**For production (best features):**
- Use 3B model
- Use full dataset
- Train for 20000-50000 steps
- Takes several hours on GPU

**Memory management:**
- 70M model: ~4GB GPU memory
- 1.3B model: ~8GB GPU memory
- 3B model: ~16GB GPU memory
- SAE training adds ~2-4GB

## Repository Structure

```
sae_genomics/
├── external/                      # Git submodules
│   ├── tahoe-x1/                 # Tahoe x1 model
│   └── SAELens/                  # SAE training framework
├── src/sae_genomics/             # Main package
│   ├── models/                   # Model adapters
│   ├── training/                 # SAE training
│   ├── validation/               # Biological validation
│   │   ├── databases/           # Database clients
│   │   ├── enrichment/          # Statistical analysis
│   │   └── visualization/       # Plots and reports
│   ├── data/                    # Data loading
│   └── utils/                   # Utilities
├── configs/                      # Configuration files
│   ├── models/                  # Model configs
│   ├── training/                # Training configs
│   └── validation/              # Validation configs
├── scripts/                      # Command-line scripts
├── notebooks/                    # Jupyter notebooks
├── tests/                       # Unit tests
└── docs/                        # Documentation
```

## Configuration

### Model Configuration

Configure tahoe x1 model size and parameters:

```yaml
# configs/models/tx1_70m.yaml
model:
  name: "tahoe-x1-70m"
  n_params: 70_000_000
  context_length: 1024
  device: "auto"
  dtype: "float32"
```

### Training Configuration

Configure SAE architecture and training:

```yaml
# configs/training/default.yaml
training:
  sae:
    d_in: 512
    d_sae: 16384
    activation: "topk"
    k: 64

  optimizer:
    lr: 3.0e-4

  batch_size: 32
  total_training_steps: 100000
```

### Validation Configuration

Configure which databases to use:

```yaml
# configs/validation/databases.yaml
databases:
  open_targets:
    enabled: true
  hpo:
    enabled: true
  # ... configure all 10 databases
```

## Biological Validation

### Feature-Disease Association

For each SAE feature:
1. Extract top N genes by activation
2. Query disease databases
3. Perform enrichment analysis (Fisher's exact test)
4. Correct for multiple testing (FDR)
5. Report significant associations

### Phenotype Enrichment

Map features to HPO terms:
- Hierarchical phenotype analysis
- Clinical feature validation
- Cross-disease phenotype patterns

### Pathway Analysis

Test features against:
- Reactome pathways
- KEGG pathways
- GO biological processes
- Gene Ontology terms

### Network Validation

Using STRING:
- Protein-protein interactions
- Functional coherence
- Network modules

## Development

### Running Tests

```bash
pytest tests/
```

### Code Quality

```bash
# Format code
black src/ tests/

# Lint
ruff check src/ tests/

# Type checking
mypy src/
```

### Pre-commit Hooks

```bash
pre-commit install
pre-commit run --all-files
```

## Roadmap

### Phase 1: Core Infrastructure (Current)
- [x] Repository structure
- [x] Configuration system
- [ ] Tahoe x1 model adapter
- [ ] SAE training pipeline
- [ ] Basic validation

### Phase 2: Validation Suite
- [ ] All 10 database clients
- [ ] Enrichment analysis
- [ ] Visualization tools
- [ ] Validation reports

### Phase 3: Production Features
- [ ] Distributed training (3B model)
- [ ] Activation caching
- [ ] Feature steering
- [ ] Advanced visualizations
- [ ] Web dashboard

### Phase 4: Research Tools
- [ ] Feature attribution
- [ ] Cross-model comparison
- [ ] Perturbation analysis
- [ ] Novel feature discovery

## Citation

If you use this repository, please cite:

```bibtex
@software{sae_genomics,
  title = {SAE Genomics: Sparse Autoencoders for Single-Cell Foundation Models},
  author = {Your Name},
  year = {2025},
  url = {https://github.com/yourusername/sae-genomics}
}
```

Also cite the underlying tools:

**Tahoe X1:**
```bibtex
@software{tahoe_x1,
  title = {Tahoe X1: Perturbation-Trained Single-Cell Foundation Model},
  author = {Tahoe Therapeutics},
  year = {2024},
  url = {https://github.com/tahoebio/tahoe-x1}
}
```

**SAELens:**
```bibtex
@software{saelens,
  title = {SAELens: Training and Analyzing Sparse Autoencoders},
  author = {Joseph Bloom and others},
  year = {2024},
  url = {https://github.com/jbloomAus/SAELens}
}
```

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Submit a pull request

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contact

For questions or issues:
- Open an issue on GitHub
- Email: your.email@example.com

## Acknowledgments

- Tahoe Therapeutics for the tahoe x1 model
- Joseph Bloom and contributors for SAELens
- Open source communities maintaining biological databases
- Research groups advancing interpretability in genomics

---

**Status:** Alpha - Active Development

Last updated: 2025-11-01
