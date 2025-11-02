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

# Lazy imports - clients will be implemented when needed
__all__ = ["OpenTargetsClient", "HPOClient", "GWASCatalogClient"]
