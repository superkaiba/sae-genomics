"""Statistical enrichment analysis for validation."""

__all__ = ["EnrichmentAnalyzer", "EnrichmentResult"]


def __getattr__(name: str):
    """Lazy imports for enrichment modules."""
    if name == "EnrichmentAnalyzer":
        from sae_genomics.validation.enrichment.fisher_exact import EnrichmentAnalyzer
        return EnrichmentAnalyzer

    if name == "EnrichmentResult":
        from sae_genomics.validation.enrichment.fisher_exact import EnrichmentResult
        return EnrichmentResult

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
