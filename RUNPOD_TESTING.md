# RunPod Testing Instructions

## Your RunPod Setup

**Instance**: 63.141.33.89:22120
**GPUs**: 3x NVIDIA H100 80GB HBM3
**Python**: 3.11.11
**OS**: Ubuntu 22.04.5

## Step-by-Step Testing

### 1. Connect to RunPod

```bash
ssh root@63.141.33.89 -p 22120
# Password: Superkaiba15!
```

### 2. Setup Repository on RunPod

Once connected, run these commands:

```bash
# Go to home directory
cd ~

# Install uv if not present
which uv || curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH=$HOME/.local/bin:$PATH

# Clone your local repository (we'll upload the git bundle)
# For now, create the structure manually
mkdir -p sae_genomics
cd sae_genomics
```

### 3. Option A: Upload via Git Bundle (From Your Mac)

**On your Mac**, run:
```bash
# Create bundle
git bundle create /tmp/sae_genomics.bundle --all

# Upload (you'll need to enter password)
scp -P 22120 /tmp/sae_genomics.bundle root@63.141.33.89:/root/
```

**On RunPod**, run:
```bash
cd ~
git clone sae_genomics.bundle sae_genomics
cd sae_genomics
git submodule update --init --recursive
```

### 3. Option B: Clone Directly (If Pushed to GitHub)

If you push to GitHub first:

**On your Mac**:
```bash
# Add remote and push
git remote add origin https://github.com/yourusername/sae-genomics.git
git push -u origin main
git push origin --tags
```

**On RunPod**:
```bash
git clone --recursive https://github.com/yourusername/sae-genomics.git
cd sae_genomics
```

### 3. Option C: Manual File Transfer (Simplest)

**On RunPod**, just create the files:

```bash
cd ~
git clone https://github.com/tahoebio/tahoe-x1.git sae_genomics/external/tahoe-x1
git clone https://github.com/jbloomAus/SAELens.git sae_genomics/external/SAELens
```

Then I'll provide you the key Python files to paste.

### 4. Install Dependencies

```bash
cd ~/sae_genomics
export PATH=$HOME/.local/bin:$PATH

# Create venv
uv venv
source .venv/bin/activate

# Install package
uv pip install scanpy
uv pip install -e .
```

### 5. Download Test Data

```bash
python scripts/download_test_data.py
```

Expected output:
```
✓ Saved to data/pbmc3k_test.h5ad
Dataset info:
  Cells: 2700
  Genes: 13714
  Cell types: 5
```

### 6. Run Pipeline Test

```bash
python scripts/run_pipeline.py \
  --data data/pbmc3k_test.h5ad \
  --model-size 70m \
  --output results/gpu_test \
  --steps 500 \
  --max-cells 500 \
  --device cuda
```

**Expected runtime**: ~5-10 minutes with H100s
**Expected output**:
```
results/gpu_test/
├── sae_final.pt
├── activations.pt
├── feature_gene_associations.json
├── feature_cell_associations.json
├── dashboard.html
└── training_metrics.json
```

### 7. Download Results (From Your Mac)

```bash
# Download the dashboard
scp -P 22120 root@63.141.33.89:~/sae_genomics/results/gpu_test/dashboard.html ./results_from_runpod.html

# Open it
open ./results_from_runpod.html
```

## Quick Copy-Paste Script for RunPod

```bash
# === RUN THIS ON RUNPOD ===

# Setup
cd ~
export PATH=$HOME/.local/bin:$PATH

# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Get submodules
mkdir -p sae_genomics/external
cd sae_genomics/external
git clone https://github.com/tahoebio/tahoe-x1.git
git clone https://github.com/jbloomAus/SAELens.git
cd ..

# Copy config and script files (I'll provide these)
# ...

# Install deps
uv venv && source .venv/bin/activate
uv pip install scanpy
uv pip install torch transformers sae-lens boto3 composer omegaconf mosaicml-streaming llm-foundry rich typer pyyaml python-dotenv pydantic

# Get test data
# (I'll provide download script)

# Run pipeline
# (I'll provide the scripts)
```

## What I Need From You

Since SSH automation is complex, the easiest path is:

**Option 1**: You manually SSH in and run the commands above

**Option 2**: Push this repo to GitHub, then clone on RunPod

**Option 3**: I create a single Python script that you can copy-paste into RunPod that sets everything up

Which would you prefer?
