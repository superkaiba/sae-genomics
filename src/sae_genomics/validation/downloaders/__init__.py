"""Database downloaders for biological validation."""

from typing import Type, Dict

__all__ = [
    "DatabaseDownloader",
    "DOWNLOADERS",
]


def __getattr__(name: str):
    """Lazy imports for downloaders."""
    if name == "DatabaseDownloader":
        from sae_genomics.validation.downloaders.base import DatabaseDownloader
        return DatabaseDownloader

    # Individual downloaders
    downloaders_map = {
        "OpenTargetsDownloader": "open_targets",
        "HPODownloader": "hpo",
        "GWASCatalogDownloader": "gwas_catalog",
        "DisGeNETDownloader": "disgenet",
        "MonarchDownloader": "monarch",
        "STRINGDownloader": "string",
        "OMIMDownloader": "omim",
        "OrphanetDownloader": "orphanet",
        "ClinVarDownloader": "clinvar",
        "GenCCDownloader": "gencc",
    }

    if name in downloaders_map:
        module_name = downloaders_map[name]
        module = __import__(
            f"sae_genomics.validation.downloaders.{module_name}",
            fromlist=[name]
        )
        return getattr(module, name)

    if name == "DOWNLOADERS":
        # Registry of all downloaders
        return {
            # Tier 1: Small, essential databases
            'gencc': __getattr__("GenCCDownloader"),
            'hpo': __getattr__("HPODownloader"),
            'string': __getattr__("STRINGDownloader"),
            # Tier 2: Medium-sized databases
            'gwas_catalog': __getattr__("GWASCatalogDownloader"),
            'orphanet': __getattr__("OrphanetDownloader"),
            'clinvar': __getattr__("ClinVarDownloader"),
            # Tier 3: Large databases
            'monarch': __getattr__("MonarchDownloader"),
            'open_targets': __getattr__("OpenTargetsDownloader"),
            # API-only databases (cache setup)
            'disgenet': __getattr__("DisGeNETDownloader"),
            'omim': __getattr__("OMIMDownloader"),
        }

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
