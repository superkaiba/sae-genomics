# Configuration Guide

SAE Genomics uses YAML configuration files for all major components.

## Configuration Files

### Model Configuration
- `configs/models/tx1_70m.yaml` - 70M parameter model (prototyping)
- `configs/models/tx1_3b.yaml` - 3B parameter model (production)

### Training Configuration
- `configs/training/default.yaml` - Default training settings

### Validation Configuration
- `configs/validation/databases.yaml` - Database access and settings

## Model Configuration

```yaml
model:
  name: "tahoe-x1-70m"
  size: "70M"

  # Architecture
  n_params: 70_000_000
  context_length: 1024
  d_model: 512
  n_layers: 12
  n_heads: 8

  # Loading
  checkpoint_path: null
  from_pretrained: true
  cache_dir: "./data/models"

  # Compute
  device: "auto"  # auto, cuda, cpu, mps
  dtype: "float32"  # float32, bfloat16, float16

  # Hook points for activation extraction
  hook_points:
    - "blocks.0.hook_mlp_out"
    - "blocks.3.hook_mlp_out"
    # ...
```

## Training Configuration

```yaml
training:
  # Experiment tracking
  experiment_name: "sae_training"
  wandb_project: "sae-genomics"

  # SAE architecture
  sae:
    d_in: 512
    d_sae: 16384
    expansion_factor: 32
    activation: "topk"
    k: 64

  # Optimization
  optimizer:
    name: "adam"
    lr: 3.0e-4
    betas: [0.9, 0.999]

  # Training loop
  batch_size: 32
  total_training_steps: 100000
  eval_every: 1000
  save_every: 5000

  # Reproducibility
  seed: 42
```

## Validation Configuration

```yaml
databases:
  # Enable/disable databases
  open_targets:
    enabled: true
    api_url: "https://api.platform.opentargets.org/api/v4/graphql"
    cache_dir: "./data/validation_cache/open_targets"
    cache_ttl: 604800  # 7 days

  # API keys (from environment)
  disgenet:
    enabled: true
    api_key: null  # Set via DISGENET_API_KEY env var

  # Local data files
  hpo:
    enabled: true
    obo_file: "./data/databases/hp.obo"
    auto_download: true

enrichment:
  # Statistical testing
  test_method: "fisher"
  multiple_testing_correction: "fdr_bh"
  significance_threshold: 0.05

  # Gene sets
  top_n_genes: 100
  min_gene_set_size: 5
  max_gene_set_size: 500
```

## Environment Variables

Set in `.env` file:

```bash
# Required
WANDB_API_KEY=your_key_here
DISGENET_API_KEY=your_key_here
OMIM_API_KEY=your_key_here

# Optional
CUDA_VISIBLE_DEVICES=0
LOG_LEVEL=INFO
```

## Command-Line Overrides

Override config values from command line:

```bash
sae-train \
    --config configs/training/default.yaml \
    --model-config configs/models/tx1_3b.yaml \
    --output ./results/experiment_1
```

## Advanced Configuration

### Distributed Training

For large models:

```yaml
model:
  distributed:
    enabled: true
    strategy: "fsdp"  # or "deepspeed"
```

### Custom Hook Points

Specify which layers to extract activations from:

```yaml
model:
  hook_points:
    - "blocks.0.hook_mlp_out"
    - "blocks.5.hook_mlp_out"
    - "blocks.10.hook_mlp_out"
```

### Validation Subsets

Validate only specific features:

```bash
sae-validate checkpoint.pt --features 0,42,100-200
```

## Configuration Best Practices

1. **Use version control** for configs
2. **Document changes** in config comments
3. **Use environment variables** for secrets
4. **Create experiment-specific configs** for reproducibility
5. **Test with small configs** before scaling up

## Next Steps

- See [Training Guide](guides/training.md) for training details
- See [Validation Guide](guides/validation.md) for validation options
- Check [API Reference](api/index.md) for programmatic configuration
