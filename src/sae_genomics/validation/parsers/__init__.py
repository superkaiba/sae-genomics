"""Database parsers for validation data."""

__all__ = [
    "DatabaseParser",
    "GeneDiseaseAssociation",
    "GenePhenotypeAssociation",
    "GenCCParser",
    "HPOParser",
    "GWASCatalogParser",
    "STRINGParser",
    "OrphanetParser",
    "ClinVarParser",
    "DisGeNETParser",
    "OMIMParser",
    "MonarchParser",
    "OpenTargetsParser",
]


def __getattr__(name: str):
    """Lazy imports for parsers."""
    if name == "DatabaseParser":
        from sae_genomics.validation.parsers.base import DatabaseParser
        return DatabaseParser

    if name == "GeneDiseaseAssociation":
        from sae_genomics.validation.parsers.base import GeneDiseaseAssociation
        return GeneDiseaseAssociation

    if name == "GenePhenotypeAssociation":
        from sae_genomics.validation.parsers.base import GenePhenotypeAssociation
        return GenePhenotypeAssociation

    # Specific parsers
    parsers_map = {
        "GenCCParser": "gencc",
        "HPOParser": "hpo",
        "GWASCatalogParser": "gwas_catalog",
        "STRINGParser": "string",
        "OrphanetParser": "orphanet",
        "ClinVarParser": "clinvar",
        "DisGeNETParser": "disgenet",
        "OMIMParser": "omim",
        "MonarchParser": "monarch",
        "OpenTargetsParser": "open_targets",
    }

    if name in parsers_map:
        module_name = parsers_map[name]
        module = __import__(
            f"sae_genomics.validation.parsers.{module_name}",
            fromlist=[name]
        )
        return getattr(module, name)

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
