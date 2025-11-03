"""Base class for database downloaders."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Optional, List
from datetime import datetime
import requests
import hashlib
import json
import gzip
import shutil

from rich.progress import (
    Progress,
    DownloadColumn,
    BarColumn,
    TimeRemainingColumn,
    TextColumn,
    TransferSpeedColumn,
)


class DatabaseDownloader(ABC):
    """Base class for database downloaders.

    Provides common functionality for downloading, verifying, and managing
    biological validation databases.
    """

    def __init__(self, output_dir: Path, **kwargs):
        """Initialize downloader.

        Args:
            output_dir: Directory to store downloaded files
            **kwargs: Additional arguments for specific downloaders
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.metadata_file = self.output_dir / 'metadata.json'
        self.kwargs = kwargs

    @abstractmethod
    def get_download_urls(self) -> List[Dict[str, str]]:
        """Get list of files to download.

        Returns:
            List of dicts with keys:
                - url: Download URL
                - path: Local file path (relative to output_dir)
                - decompress: Whether to decompress after download (default False)
                - verify_checksum: Optional checksum to verify (default None)
        """
        pass

    @abstractmethod
    def get_latest_version(self) -> str:
        """Get the latest version identifier for this database.

        Returns:
            Version string (e.g., "2025-10-22", "v12.0", etc.)
        """
        pass

    def download_file(
        self,
        url: str,
        output_path: Path,
        chunk_size: int = 8192,
        show_progress: bool = True,
    ) -> Path:
        """Download a file with progress bar.

        Args:
            url: URL to download
            output_path: Local path to save file
            chunk_size: Download chunk size in bytes
            show_progress: Whether to show progress bar

        Returns:
            Path to downloaded file

        Raises:
            requests.HTTPError: If download fails
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Start download
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))

        if show_progress and total_size > 0:
            with Progress(
                TextColumn("[bold blue]{task.description}"),
                BarColumn(),
                "[progress.percentage]{task.percentage:>3.0f}%",
                DownloadColumn(),
                TransferSpeedColumn(),
                TimeRemainingColumn(),
            ) as progress:
                task = progress.add_task(
                    f"Downloading {output_path.name}",
                    total=total_size
                )

                with open(output_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            progress.update(task, advance=len(chunk))
        else:
            # No progress bar for small files or unknown size
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)

        return output_path

    def decompress_file(
        self,
        gz_path: Path,
        output_path: Optional[Path] = None
    ) -> Path:
        """Decompress a gzip file.

        Args:
            gz_path: Path to gzipped file
            output_path: Output path (defaults to gz_path without .gz extension)

        Returns:
            Path to decompressed file
        """
        if output_path is None:
            # Remove .gz extension
            if gz_path.suffix == '.gz':
                output_path = gz_path.with_suffix('')
            else:
                output_path = gz_path.parent / f"{gz_path.name}.decompressed"

        print(f"Decompressing {gz_path.name}...")
        with gzip.open(gz_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        return output_path

    def compute_checksum(
        self,
        file_path: Path,
        algorithm: str = 'sha256'
    ) -> str:
        """Compute file checksum.

        Args:
            file_path: Path to file
            algorithm: Hash algorithm (sha256, md5, etc.)

        Returns:
            Hexadecimal checksum string
        """
        hasher = hashlib.new(algorithm)
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(65536), b''):
                hasher.update(chunk)
        return hasher.hexdigest()

    def download(self) -> Dict:
        """Download all files for this database.

        Returns:
            Metadata dict with version, download date, and file info

        Raises:
            Exception: If download or verification fails
        """
        files = self.get_download_urls()
        metadata = {
            'version': self.get_latest_version(),
            'download_date': datetime.now().isoformat(),
            'database_name': self.__class__.__name__.replace('Downloader', '').lower(),
            'files': []
        }

        print(f"\nDownloading {metadata['database_name']} version {metadata['version']}")
        print(f"Output directory: {self.output_dir}")
        print(f"Files to download: {len(files)}\n")

        for file_info in files:
            url = file_info['url']
            local_path = self.output_dir / file_info['path']

            # Download
            try:
                downloaded_path = self.download_file(url, local_path)
            except Exception as e:
                print(f"✗ Failed to download {file_info['path']}: {e}")
                raise

            # Decompress if needed
            if file_info.get('decompress', False):
                try:
                    decompressed_path = self.decompress_file(downloaded_path)
                    final_path = decompressed_path
                except Exception as e:
                    print(f"✗ Failed to decompress {downloaded_path.name}: {e}")
                    raise
            else:
                final_path = downloaded_path

            # Compute checksum
            checksum = self.compute_checksum(final_path)

            # Verify checksum if provided
            expected_checksum = file_info.get('verify_checksum')
            if expected_checksum and checksum != expected_checksum:
                raise ValueError(
                    f"Checksum mismatch for {final_path.name}\n"
                    f"Expected: {expected_checksum}\n"
                    f"Got: {checksum}"
                )

            size_mb = final_path.stat().st_size / (1024 * 1024)

            # Store the actual file path (decompressed if applicable)
            metadata['files'].append({
                'name': str(final_path.relative_to(self.output_dir)),
                'size_mb': round(size_mb, 2),
                'checksum': checksum,
                'decompressed': file_info.get('decompress', False)
            })

            print(f"✓ {final_path.name} ({size_mb:.1f} MB)")

        # Database-specific post-processing
        self.post_process()

        # Save metadata
        self.save_metadata(metadata)

        print(f"\n✓ Download complete for {metadata['database_name']}")
        return metadata

    def post_process(self):
        """Optional post-processing after download.

        Override in subclasses if needed (e.g., to build indexes,
        extract archives, etc.).
        """
        pass

    def verify(self) -> bool:
        """Verify that all downloaded files are present and valid.

        Returns:
            True if all files are valid, False otherwise
        """
        if not self.metadata_file.exists():
            return False

        try:
            with open(self.metadata_file) as f:
                metadata = json.load(f)
        except Exception:
            return False

        # Check all files exist and have correct checksums
        for file_info in metadata.get('files', []):
            file_path = self.output_dir / file_info['name']
            if not file_path.exists():
                return False

            # Verify checksum
            checksum = self.compute_checksum(file_path)
            if checksum != file_info['checksum']:
                return False

        return True

    def is_downloaded(self) -> bool:
        """Check if database is already downloaded and verified.

        Returns:
            True if database is downloaded and valid
        """
        return self.metadata_file.exists() and self.verify()

    def get_metadata(self) -> Dict:
        """Load metadata from file.

        Returns:
            Metadata dict

        Raises:
            FileNotFoundError: If metadata file doesn't exist
        """
        if not self.metadata_file.exists():
            raise FileNotFoundError(f"Metadata file not found: {self.metadata_file}")

        with open(self.metadata_file) as f:
            return json.load(f)

    def save_metadata(self, metadata: Dict):
        """Save metadata to file.

        Args:
            metadata: Metadata dict to save
        """
        with open(self.metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        # Also create a version.txt file for easy checking
        version_file = self.output_dir / 'version.txt'
        with open(version_file, 'w') as f:
            f.write(f"{metadata['version']}\n")
            f.write(f"Downloaded: {metadata['download_date']}\n")
