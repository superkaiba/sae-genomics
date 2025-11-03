"""STRING protein-protein interaction database downloader."""

from pathlib import Path
from typing import List, Dict

from .base import DatabaseDownloader


class STRINGDownloader(DatabaseDownloader):
    """Download STRING protein-protein interaction data for human.

    STRING provides comprehensive protein-protein interaction networks
    with confidence scores from multiple evidence channels.

    Database info:
        - Size: ~80-150 MB compressed, ~600 MB-1 GB decompressed
        - Format: TSV (tab-separated, gzipped)
        - License: CC BY 4.0
        - Updates: Annual major releases
        - Current version: 12.0
        - URL: https://string-db.org
    """

    BASE_URL = "https://stringdb-downloads.org/download"
    VERSION = "12.0"
    SPECIES = "9606"  # Human NCBI taxonomy ID

    def __init__(self, output_dir: Path, detailed: bool = True, **kwargs):
        """Initialize STRING downloader.

        Args:
            output_dir: Directory to store downloaded files
            detailed: If True, download detailed scores by channel (larger file)
            **kwargs: Additional arguments
        """
        super().__init__(output_dir, **kwargs)
        self.detailed = detailed

    def get_latest_version(self) -> str:
        """Get STRING version.

        Returns:
            Version string (e.g., "12.0")
        """
        return self.VERSION

    def get_download_urls(self) -> List[Dict[str, str]]:
        """Get STRING download URLs.

        Returns:
            List of files to download
        """
        files = []

        # Protein interaction links
        if self.detailed:
            # Detailed file with scores by channel (experimental, database, etc.)
            filename = f"{self.SPECIES}.protein.links.detailed.v{self.VERSION}.txt.gz"
        else:
            # Basic file with combined scores only
            filename = f"{self.SPECIES}.protein.links.v{self.VERSION}.txt.gz"

        files.append({
            'url': f"{self.BASE_URL}/protein.links.detailed.v{self.VERSION}/{filename}",
            'path': filename,
            'decompress': True,  # Decompress for faster queries
        })

        # Protein information (names, annotations)
        info_filename = f"{self.SPECIES}.protein.info.v{self.VERSION}.txt.gz"
        files.append({
            'url': f"{self.BASE_URL}/protein.info.v{self.VERSION}/{info_filename}",
            'path': info_filename,
            'decompress': True,
        })

        return files
