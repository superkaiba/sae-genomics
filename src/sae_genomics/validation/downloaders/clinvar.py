"""ClinVar clinical variant database downloader."""

from pathlib import Path
from typing import List, Dict
from datetime import datetime
import requests

from .base import DatabaseDownloader


class ClinVarDownloader(DatabaseDownloader):
    """Download ClinVar clinical variant interpretations.

    ClinVar provides clinical significance annotations for genetic variants,
    aggregated from clinical testing labs and research studies.

    Database info:
        - Size: ~200-400 MB compressed, ~1.5-2 GB decompressed
        - Format: VCF (Variant Call Format) or TSV
        - License: Public domain (US government data)
        - Updates: Monthly (first Thursday) + weekly XMLs
        - URL: https://www.ncbi.nlm.nih.gov/clinvar/
    """

    FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar"

    def __init__(self, output_dir: Path, format: str = 'tsv', genome_build: str = 'GRCh38', **kwargs):
        """Initialize ClinVar downloader.

        Args:
            output_dir: Directory to store downloaded files
            format: 'vcf' for VCF format or 'tsv' for variant summary table (default: 'tsv')
            genome_build: 'GRCh38' or 'GRCh37' (default: 'GRCh38')
            **kwargs: Additional arguments
        """
        super().__init__(output_dir, **kwargs)
        self.format = format
        self.genome_build = genome_build

    def get_latest_version(self) -> str:
        """Get latest version (uses download date).

        ClinVar doesn't expose version numbers easily, so we use download date.

        Returns:
            Current date as version string
        """
        return datetime.now().strftime('%Y-%m-%d')

    def get_download_urls(self) -> List[Dict[str, str]]:
        """Get ClinVar download URLs.

        Returns:
            List of files to download
        """
        if self.format == 'vcf':
            # VCF format (more complete but harder to parse)
            return [
                {
                    'url': f"{self.FTP_BASE}/vcf_{self.genome_build}/clinvar.vcf.gz",
                    'path': f'clinvar_{self.genome_build}.vcf.gz',
                    'decompress': True,
                },
                # Index file
                {
                    'url': f"{self.FTP_BASE}/vcf_{self.genome_build}/clinvar.vcf.gz.tbi",
                    'path': f'clinvar_{self.genome_build}.vcf.gz.tbi',
                    'decompress': False,
                }
            ]
        else:
            # TSV format (easier for gene-level queries)
            return [
                {
                    'url': f"{self.FTP_BASE}/tab_delimited/variant_summary.txt.gz",
                    'path': 'variant_summary.txt.gz',
                    'decompress': True,
                }
            ]
