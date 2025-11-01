# Notebooks

Interactive Jupyter notebooks for exploring SAE Genomics.

## Exploratory

### 01_data_exploration.ipynb
- Load and explore tahoe x1 model
- Visualize single-cell data
- Understand model activations

### 02_sae_training.ipynb
- Train SAEs interactively
- Monitor training metrics
- Visualize learned features

### 03_feature_analysis.ipynb
- Analyze SAE features
- Identify top activated genes
- Explore feature relationships

## Validation

### 04_validation_results.ipynb
- Load validation results
- Visualize enrichment analysis
- Generate custom plots

### 05_database_exploration.ipynb
- Query biological databases
- Test enrichment methods
- Compare database results

## Usage

Start Jupyter Lab:
```bash
jupyter lab
```

Or use VS Code with Jupyter extension.

## Requirements

Install notebook dependencies:
```bash
uv pip install -e ".[notebooks]"
```
