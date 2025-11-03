"""GenCC (Gene Curation Coalition) database downloader."""

from pathlib import Path
from typing import List, Dict
from datetime import datetime

from .base import DatabaseDownloader


class GenCCDownloader(DatabaseDownloader):
    """Download GenCC gene-disease validity classifications.

    GenCC provides curated gene-disease associations with validity
    classifications from multiple expert panels (ClinGen, Genomics England, etc.).

    Database info:
        - Size: ~5 MB
        - Format: CSV
        - License: CC0 (public domain)
        - Updates: Continuous
        - URL: https://search.thegencc.org/download
    """

    BASE_URL = "https://search.thegencc.org/download/action/submissions-export-csv"

    def get_latest_version(self) -> str:
        """Get latest version (uses download date as version).

        Returns:
            Current date as version string
        """
        return datetime.now().strftime('%Y-%m-%d')

    def get_download_urls(self) -> List[Dict[str, str]]:
        """Get GenCC CSV download URL.

        Returns:
            List with single CSV file to download
        """
        return [
            {
                'url': self.BASE_URL,
                'path': 'gencc-submissions.csv',
                'decompress': False,
            }
        ]
