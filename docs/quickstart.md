# Quick Start Guide

This guide will walk you through training your first SAE and validating its features.

## 1. Training an SAE

### Using the Command Line

```bash
sae-train --config configs/training/default.yaml --output ./results/my_first_sae
```

### Using Python

```python
from sae_genomics.training import SAETrainer
from sae_genomics.models import TahoeModelAdapter

# Load tahoe x1 model (start with 70M for prototyping)
model = TahoeModelAdapter.from_config("configs/models/tx1_70m.yaml")

# Initialize trainer
trainer = SAETrainer(config="configs/training/default.yaml")

# Train the SAE
trainer.train(model, output_dir="./results/my_first_sae")
```

### In a Notebook

See `notebooks/exploratory/02_sae_training.ipynb` for an interactive example.

## 2. Validating Features

Once training is complete, validate the learned features:

```bash
sae-validate results/my_first_sae/checkpoints/final.pt \
    --config configs/validation/databases.yaml \
    --output results/my_first_sae/validation
```

This will:
1. Load the trained SAE
2. Extract feature activations
3. Query all 10 biological databases
4. Perform enrichment analysis
5. Generate validation report

## 3. Analyzing Results

### View Validation Report

```bash
open results/my_first_sae/validation/report.html
```

### Explore in Python

```python
from sae_genomics.validation import load_validation_results

results = load_validation_results("results/my_first_sae/validation")

# Get top disease associations for feature 42
feature_42 = results.get_feature(42)
print(feature_42.disease_enrichment.head())

# Plot enrichment
feature_42.plot_enrichment(save_path="feature_42_enrichment.png")
```

### Interactive Notebook

See `notebooks/validation/04_validation_results.ipynb` for detailed analysis.

## 4. Customization

### Adjust Model Size

For production, use the 3B model:

```yaml
# configs/training/production.yaml
training:
  model_config: "configs/models/tx1_3b.yaml"
  # ... adjust batch size, etc.
```

### Configure SAE Architecture

```yaml
# configs/training/default.yaml
sae:
  d_in: 512
  d_sae: 16384  # Dictionary size
  activation: "topk"
  k: 64  # Number of active features
```

### Select Databases

```yaml
# configs/validation/databases.yaml
databases:
  open_targets:
    enabled: true
  hpo:
    enabled: true
  gwas_catalog:
    enabled: true
  # Disable others if needed
  disgenet:
    enabled: false
```

## Next Steps

- Read the [Training Guide](guides/training.md) for advanced options
- Learn about [Biological Validation](guides/validation.md)
- Explore [Configuration Options](configuration.md)
- Check out [Example Notebooks](../notebooks/README.md)
