"""Biological validation tools for SAE features."""

# Lazy imports - modules will be loaded when accessed
__all__ = ["databases", "enrichment", "visualization", "parsers", "LocalValidationClient"]


def __getattr__(name: str):
    """Lazy imports for validation modules."""
    if name == "LocalValidationClient":
        from sae_genomics.validation.local_client import LocalValidationClient
        return LocalValidationClient

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
