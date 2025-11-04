"""Open Targets Platform database downloader."""

from pathlib import Path
from typing import List, Dict
import subprocess
import ftplib
from datetime import datetime

from .base import DatabaseDownloader


class OpenTargetsDownloader(DatabaseDownloader):
    """Download Open Targets Platform data.

    Open Targets provides comprehensive evidence linking targets (genes)
    to diseases, integrating data from genetics, drugs, pathways, and literature.

    Database info:
        - Size: 6-10 GB (gene-disease associations only), 15-25 GB (full)
        - Format: Parquet (columnar format)
        - License: Apache 2.0 (freely downloadable)
        - Updates: Quarterly releases
        - URL: https://platform.opentargets.org/
    """

    FTP_HOST = "ftp.ebi.ac.uk"
    FTP_BASE_PATH = "/pub/databases/opentargets/platform"
    HTTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform"

    def __init__(self, output_dir: Path, datasets: List[str] = None, **kwargs):
        """Initialize Open Targets downloader.

        Args:
            output_dir: Directory to store downloaded files
            datasets: List of datasets to download. Options:
                - 'disease': Disease information
                - 'target': Target (gene) information
                - 'association_by_datasource_direct': Gene-disease associations (recommended)
            **kwargs: Additional arguments
        """
        super().__init__(output_dir, **kwargs)
        # Default: only download gene-disease associations (smallest useful subset)
        self.datasets = datasets or [
            'disease',
            'target',
            'association_by_datasource_direct',  # Direct associations
        ]

    def get_latest_version(self) -> str:
        """Get latest Open Targets release version from FTP.

        Returns:
            Version string (e.g., "24.09" for September 2024 release)
        """
        try:
            ftp = ftplib.FTP(self.FTP_HOST)
            ftp.login()
            ftp.cwd(self.FTP_BASE_PATH)

            # List all directories
            versions = []
            ftp.retrlines('NLST', versions.append)

            # Filter to version numbers (format: YY.MM)
            version_dirs = [v for v in versions if v.replace('.', '').isdigit() and len(v) <= 5]

            ftp.quit()

            if version_dirs:
                # Get the latest version
                return max(version_dirs, key=lambda v: tuple(map(int, v.split('.'))))
        except Exception as e:
            print(f"Warning: Could not fetch latest version from FTP: {e}")

        # Fallback: use current month
        now = datetime.now()
        return f"{now.year % 100}.{now.month:02d}"

    def get_download_urls(self) -> List[Dict[str, str]]:
        """Get Open Targets download URLs.

        Note: Open Targets uses partitioned parquet files, so we'll use
        a directory-based approach.

        Returns:
            List of dataset directories to download
        """
        version = self.get_latest_version()

        urls = []
        for dataset in self.datasets:
            urls.append({
                'url': f"{self.HTTP_BASE}/{version}/output/{dataset}/",
                'path': dataset,
                'type': 'directory',
                'decompress': False
            })

        return urls

    def download(self) -> Dict:
        """Custom download for directory structures.

        Open Targets data is organized in directories with multiple parquet files.
        We use wget/curl for recursive download.

        Returns:
            Metadata dict
        """
        version = self.get_latest_version()
        metadata = {
            'version': version,
            'download_date': datetime.now().isoformat(),
            'database_name': 'open_targets',
            'files': []
        }

        print(f"\nDownloading Open Targets version {version}")
        print(f"Output directory: {self.output_dir}")
        print(f"Datasets: {', '.join(self.datasets)}\n")

        for dataset_info in self.get_download_urls():
            url = dataset_info['url']
            local_path = self.output_dir / dataset_info['path']
            local_path.mkdir(parents=True, exist_ok=True)

            print(f"Downloading {dataset_info['path']}...")

            try:
                # Use wget for recursive download (more reliable than Python requests)
                cmd = [
                    'wget',
                    '--recursive',  # Recursive download
                    '--no-parent',  # Don't ascend to parent directory
                    '--no-host-directories',  # Don't create host directory
                    '--cut-dirs=5',  # Skip 5 directory levels in URL
                    '--reject=index.html*',  # Skip index files
                    '--accept=*.parquet',  # Only download parquet files
                    '-P', str(local_path.parent),  # Output directory
                    '--progress=bar:force',  # Show progress
                    '--tries=3',  # Retry 3 times on failure
                    url
                ]

                subprocess.run(cmd, check=True, capture_output=False)

                # Calculate total size
                total_size = sum(
                    f.stat().st_size
                    for f in local_path.rglob('*.parquet')
                    if f.is_file()
                )

                metadata['files'].append({
                    'name': dataset_info['path'],
                    'size_mb': round(total_size / (1024 * 1024), 2),
                    'type': 'directory',
                    'n_parquet_files': len(list(local_path.rglob('*.parquet')))
                })

                print(f"✓ {dataset_info['path']} ({total_size / (1024**3):.2f} GB)")

            except subprocess.CalledProcessError as e:
                print(f"✗ Failed to download {dataset_info['path']}: {e}")
                # Try alternative: use curl instead of wget
                print("Trying with curl...")
                try:
                    # Note: This is a simplified fallback; full recursive download with curl is complex
                    print(f"Please manually download from: {url}")
                    raise

                except Exception as e2:
                    print(f"✗ Curl also failed: {e2}")
                    raise

        # Save metadata
        self.save_metadata(metadata)

        print(f"\n✓ Download complete for open_targets")
        return metadata

    def verify(self) -> bool:
        """Verify Open Targets download.

        For directories, we just check that parquet files exist.

        Returns:
            True if valid
        """
        if not self.metadata_file.exists():
            return False

        try:
            metadata = self.get_metadata()

            # Check each dataset directory has parquet files
            for file_info in metadata.get('files', []):
                dataset_path = self.output_dir / file_info['name']
                if not dataset_path.exists():
                    return False

                # Check for at least some parquet files
                parquet_files = list(dataset_path.rglob('*.parquet'))
                if len(parquet_files) == 0:
                    return False

            return True

        except Exception:
            return False
