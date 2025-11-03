"""Orphanet rare disease database downloader."""

from pathlib import Path
from typing import List, Dict
from datetime import datetime

from .base import DatabaseDownloader


class OrphanetDownloader(DatabaseDownloader):
    """Download Orphanet rare disease data.

    Orphanet provides comprehensive information on rare diseases and
    orphan drugs, with gene-disease associations.

    Database info:
        - Size: ~50 MB (all XML files)
        - Format: XML
        - License: CC BY 4.0
        - Updates: Biannually (July, December)
        - URL: http://www.orphadata.org/
    """

    BASE_URL = "http://www.orphadata.com/data/xml"

    # Key data products
    PRODUCTS = {
        'en_product1.xml': 'Rare diseases inventory',
        'en_product4.xml': 'Disease-HPO associations',
        'en_product6.xml': 'Disease-gene associations (KEY FILE)',
        'en_product9_ages.xml': 'Disease age of onset',
    }

    def get_latest_version(self) -> str:
        """Get latest version (uses download date).

        Orphanet releases twice yearly, but we use download date
        since there's no version endpoint.

        Returns:
            Current date as version string
        """
        return datetime.now().strftime('%Y-%m-%d')

    def get_download_urls(self) -> List[Dict[str, str]]:
        """Get Orphanet XML file URLs.

        Returns:
            List of XML files to download
        """
        return [
            {
                'url': f"{self.BASE_URL}/{filename}",
                'path': filename,
                'decompress': False,
            }
            for filename in self.PRODUCTS.keys()
        ]
