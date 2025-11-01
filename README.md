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

### 1. Configure Your Experiment

Choose a model configuration:
```yaml
# configs/models/tx1_70m.yaml - for prototyping
# configs/models/tx1_3b.yaml - for production
```

Customize training:
```yaml
# configs/training/default.yaml
```

### 2. Train an SAE

**Using CLI:**
```bash
sae-train --config configs/training/default.yaml --output ./results/my_experiment
```

**Using script:**
```bash
python scripts/train_sae.py --config configs/training/default.yaml
```

**In Python:**
```python
from sae_genomics.training import SAETrainer
from sae_genomics.models import TahoeModelAdapter

# Load model
model = TahoeModelAdapter.from_config("configs/models/tx1_70m.yaml")

# Initialize trainer
trainer = SAETrainer(config="configs/training/default.yaml")

# Train
trainer.train(model)
```

### 3. Validate Features

**Using CLI:**
```bash
sae-validate results/checkpoints/sae_final.pt --config configs/validation/databases.yaml
```

**Using script:**
```bash
python scripts/validate_features.py results/checkpoints/sae_final.pt
```

**In Python:**
```python
from sae_genomics.validation import FeatureValidator

validator = FeatureValidator(config="configs/validation/databases.yaml")
results = validator.validate("results/checkpoints/sae_final.pt")
validator.save_report(results, "results/validation/report.html")
```

### 4. Explore in Notebooks

```bash
jupyter lab notebooks/exploratory/
```

Example notebooks:
- `01_data_exploration.ipynb` - Explore tahoe x1 activations
- `02_sae_training.ipynb` - Interactive SAE training
- `03_feature_analysis.ipynb` - Analyze learned features
- `04_validation_results.ipynb` - Visualize validation results

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
