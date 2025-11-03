"""Monarch Initiative knowledge graph downloader."""

from pathlib import Path
from typing import List, Dict
from datetime import datetime
import requests

from .base import DatabaseDownloader


class MonarchDownloader(DatabaseDownloader):
    """Download Monarch Initiative knowledge graph.

    Monarch integrates data across multiple sources (HPO, Orphanet, OMIM, etc.)
    to provide cross-species phenotype and disease associations.

    Database info:
        - Size: ~800 MB compressed, ~2-3 GB decompressed
        - Format: DuckDB (recommended), TSV, SQLite, Neo4j, RDF
        - License: CC0 (public domain)
        - Updates: Monthly
        - URL: https://monarchinitiative.org/
    """

    BASE_URL = "https://data.monarchinitiative.org/monarch-kg"

    def __init__(self, output_dir: Path, format: str = 'duckdb', **kwargs):
        """Initialize Monarch downloader.

        Args:
            output_dir: Directory to store downloaded files
            format: 'duckdb' (recommended), 'sqlite', 'tsv', or 'neo4j'
            **kwargs: Additional arguments
        """
        super().__init__(output_dir, **kwargs)
        self.format = format

    def get_latest_version(self) -> str:
        """Get latest Monarch KG version.

        Returns:
            Version string (e.g., "2025-10" for October 2025 release)
        """
        # Try to get latest release from directory listing
        try:
            # The latest/ symlink points to the most recent release
            response = requests.head(f"{self.BASE_URL}/latest/", timeout=10, allow_redirects=True)
            # Extract version from redirect URL
            final_url = response.url
            if 'monarch-kg-' in final_url:
                version = final_url.split('monarch-kg-')[-1].rstrip('/')
                return version
        except Exception:
            pass

        # Fallback: use current month
        return datetime.now().strftime('%Y-%m')

    def get_download_urls(self) -> List[Dict[str, str]]:
        """Get Monarch KG download URLs.

        Returns:
            List of files to download
        """
        if self.format == 'duckdb':
            # DuckDB format (fastest for queries)
            return [
                {
                    'url': f"{self.BASE_URL}/latest/monarch-kg.duckdb.gz",
                    'path': 'monarch-kg.duckdb.gz',
                    'decompress': True,
                }
            ]
        elif self.format == 'sqlite':
            # SQLite format
            return [
                {
                    'url': f"{self.BASE_URL}/latest/monarch-kg.db.gz",
                    'path': 'monarch-kg.db.gz',
                    'decompress': True,
                }
            ]
        elif self.format == 'tsv':
            # TSV KGX format (nodes and edges separately)
            return [
                {
                    'url': f"{self.BASE_URL}/latest/monarch-kg-denormalized-edges.tsv.gz",
                    'path': 'monarch-kg-edges.tsv.gz',
                    'decompress': True,
                },
                {
                    'url': f"{self.BASE_URL}/latest/monarch-kg-denormalized-nodes.tsv.gz",
                    'path': 'monarch-kg-nodes.tsv.gz',
                    'decompress': True,
                }
            ]
        else:
            raise ValueError(f"Unsupported format: {self.format}. Use 'duckdb', 'sqlite', or 'tsv'")
