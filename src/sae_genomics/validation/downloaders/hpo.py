"""Human Phenotype Ontology (HPO) database downloader."""

from pathlib import Path
from typing import List, Dict
import requests

from .base import DatabaseDownloader


class HPODownloader(DatabaseDownloader):
    """Download Human Phenotype Ontology data.

    HPO provides a standardized vocabulary of phenotypic abnormalities
    encountered in human disease, with gene-phenotype associations.

    Database info:
        - Size: ~75 MB
        - Format: OBO (ontology), TSV (annotations)
        - License: CC0 (public domain)
        - Updates: Monthly (around the 22nd)
        - URL: https://github.com/obophenotype/human-phenotype-ontology
    """

    API_URL = "https://api.github.com/repos/obophenotype/human-phenotype-ontology/releases/latest"
    BASE_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download"

    # Files to download
    FILES = [
        'hp-base.obo',              # Core ontology (OBO format)
        'hp-base.json',             # Core ontology (JSON format)
        'genes_to_phenotype.txt',   # Gene-to-phenotype associations
        'genes_to_disease.txt',     # Gene-to-disease via phenotypes
        'phenotype_to_genes.txt',   # Reverse mapping
    ]

    def get_latest_version(self) -> str:
        """Get latest HPO release version from GitHub.

        Returns:
            Version tag (e.g., "v2025-10-22")

        Raises:
            requests.HTTPError: If API request fails
        """
        response = requests.get(self.API_URL, timeout=10)
        response.raise_for_status()
        return response.json()['tag_name']

    def get_download_urls(self) -> List[Dict[str, str]]:
        """Get list of HPO files to download.

        Returns:
            List of file specifications
        """
        return [
            {
                'url': f"{self.BASE_URL}/{filename}",
                'path': filename,
                'decompress': False
            }
            for filename in self.FILES
        ]
