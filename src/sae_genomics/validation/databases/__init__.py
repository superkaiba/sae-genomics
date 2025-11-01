"""Database clients for biological validation.

Includes clients for:
- Open Targets Platform
- Human Phenotype Ontology (HPO)
- GWAS Catalog
- DisGeNET
- Monarch Initiative
- STRING
- OMIM
- Orphanet
- ClinVar
- GenCC
"""

from sae_genomics.validation.databases.open_targets import OpenTargetsClient
from sae_genomics.validation.databases.hpo import HPOClient
from sae_genomics.validation.databases.gwas import GWASCatalogClient

__all__ = ["OpenTargetsClient", "HPOClient", "GWASCatalogClient"]
