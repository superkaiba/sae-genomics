"""SAE training utilities and pipelines."""

# Lazy imports to avoid early dependency loading
__all__ = ["SAETrainer", "TrainingConfig"]

def __getattr__(name):
    if name == "SAETrainer":
        from sae_genomics.training.trainer import SAETrainer
        return SAETrainer
    elif name == "TrainingConfig":
        from sae_genomics.training.config import TrainingConfig
        return TrainingConfig
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
