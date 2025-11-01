# Installation Guide

## Prerequisites

- Python 3.10 or higher
- [uv](https://github.com/astral-sh/uv) package manager
- Git with submodule support
- (Optional) CUDA-capable GPU for training large models

## Installation Steps

### 1. Install uv

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### 2. Clone the Repository

```bash
git clone --recursive https://github.com/yourusername/sae-genomics.git
cd sae-genomics
```

If you already cloned without `--recursive`:
```bash
git submodule update --init --recursive
```

### 3. Create Virtual Environment

```bash
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

### 4. Install Dependencies

**For production use:**
```bash
uv pip install -e .
```

**For development:**
```bash
uv pip install -e ".[dev]"
```

**For notebook usage:**
```bash
uv pip install -e ".[notebooks]"
```

**For full installation (all features):**
```bash
uv pip install -e ".[all]"
```

### 5. Set Up Environment Variables

```bash
cp .env.example .env
```

Edit `.env` and add your API keys:
- `WANDB_API_KEY` - Get from [wandb.ai](https://wandb.ai)
- `DISGENET_API_KEY` - Register at [DisGeNET](https://www.disgenet.org/)
- `OMIM_API_KEY` - Apply at [OMIM](https://www.omim.org/api)

### 6. Verify Installation

```bash
python -c "import sae_genomics; print(sae_genomics.__version__)"
sae-train --help
```

## Troubleshooting

### Submodule Issues

If submodules didn't clone correctly:
```bash
git submodule update --init --recursive --force
```

### CUDA Issues

Ensure PyTorch is installed with CUDA support:
```bash
uv pip install torch --index-url https://download.pytorch.org/whl/cu121
```

### Permission Issues

On Linux/Mac, you may need to make scripts executable:
```bash
chmod +x scripts/*.py
```

## Next Steps

- Read the [Quick Start Guide](quickstart.md)
- Review [Configuration Options](configuration.md)
- Explore [Example Notebooks](../notebooks/README.md)
