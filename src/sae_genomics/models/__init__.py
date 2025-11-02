"""Model loading and adaptation utilities for tahoe x1."""

# Lazy imports to avoid early dependency loading
__all__ = ["TahoeModelAdapter", "ModelConfig"]

def __getattr__(name):
    if name == "TahoeModelAdapter":
        from sae_genomics.models.tahoe_adapter import TahoeModelAdapter
        return TahoeModelAdapter
    elif name == "ModelConfig":
        from sae_genomics.models.model_config import ModelConfig
        return ModelConfig
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
