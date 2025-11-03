"""SAE Genomics: Sparse Autoencoder training for tahoe x1 single-cell foundation model."""

__version__ = "0.1.0"

# Use lazy imports to avoid circular dependencies
__all__ = ["models", "training", "validation", "utils", "__version__"]


def __getattr__(name):
    """Lazy load submodules to avoid circular imports."""
    if name in ["models", "training", "validation", "utils"]:
        import importlib
        module = importlib.import_module(f"sae_genomics.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
