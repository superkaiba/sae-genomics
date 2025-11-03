"""GWAS Catalog database downloader."""

from pathlib import Path
from typing import List, Dict
from datetime import datetime

from .base import DatabaseDownloader


class GWASCatalogDownloader(DatabaseDownloader):
    """Download GWAS Catalog gene-trait associations.

    The GWAS Catalog provides a curated collection of all published
    genome-wide association studies with SNP-trait associations.

    Database info:
        - Size: ~200-400 MB
        - Format: TSV (tab-separated values)
        - License: EMBL-EBI Terms of Use (free for academic)
        - Updates: Weekly (typically Sunday/Monday)
        - URL: https://www.ebi.ac.uk/gwas/
    """

    # Use FTP for reliable access
    FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest"

    def get_latest_version(self) -> str:
        """Get latest version (uses download date as version).

        Returns:
            Current date as version string
        """
        return datetime.now().strftime('%Y-%m-%d')

    def get_download_urls(self) -> List[Dict[str, str]]:
        """Get GWAS Catalog download URLs.

        Returns:
            List of file specifications
        """
        return [
            {
                'url': f"{self.FTP_BASE}/gwas-catalog-associations_ontology-annotated.tsv",
                'path': 'gwas-catalog-associations.tsv',
                'decompress': False,
            },
        ]
