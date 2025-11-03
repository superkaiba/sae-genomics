"""DisGeNET API cache setup downloader."""

from pathlib import Path
from typing import List, Dict
from datetime import datetime
import os

from .base import DatabaseDownloader


class DisGeNETDownloader(DatabaseDownloader):
    """Set up DisGeNET API cache structure.

    DisGeNET does not allow bulk downloads for free academic users.
    This "downloader" creates the cache structure and documentation
    for API-based access with local caching.

    Database info:
        - Size: API cache (grows with usage)
        - Format: JSON (cached API responses)
        - License: Academic license required (free but restricted)
        - Updates: Real-time via API
        - URL: https://www.disgenet.org/
    """

    API_BASE = "https://www.disgenet.org/api"
    SIGNUP_URL = "https://www.disgenet.org/academic-apply"

    def __init__(self, output_dir: Path, api_key: str = None, **kwargs):
        """Initialize DisGeNET cache setup.

        Args:
            output_dir: Directory for cache structure
            api_key: DisGeNET API key (optional, can be set later)
            **kwargs: Additional arguments
        """
        super().__init__(output_dir, **kwargs)
        self.api_key = api_key or os.getenv('DISGENET_API_KEY')

    def get_latest_version(self) -> str:
        """Get version (API-based, returns current date).

        Returns:
            Current date as version
        """
        return datetime.now().strftime('%Y-%m-%d')

    def get_download_urls(self) -> List[Dict[str, str]]:
        """No downloads for API-only database.

        Returns:
            Empty list
        """
        return []

    def download(self) -> Dict:
        """Create cache directory structure and documentation.

        Returns:
            Metadata dict
        """
        print(f"\nSetting up DisGeNET API cache...")
        print(f"Output directory: {self.output_dir}")

        # Create cache directories
        cache_dir = self.output_dir / 'cache'
        cache_dir.mkdir(parents=True, exist_ok=True)

        (cache_dir / 'gene_queries').mkdir(exist_ok=True)
        (cache_dir / 'disease_queries').mkdir(exist_ok=True)
        (cache_dir / 'variant_queries').mkdir(exist_ok=True)

        # Save API key if provided
        if self.api_key:
            with open(self.output_dir / 'api_key.txt', 'w') as f:
                f.write(self.api_key)
            print("✓ API key saved")
        else:
            print("⚠ No API key provided. You'll need to obtain one from:")
            print(f"  {self.SIGNUP_URL}")

        # Create README with instructions
        readme_content = f"""# DisGeNET API Cache

## About DisGeNET

DisGeNET integrates gene-disease associations from various sources including
GWAS, animal models, and literature mining.

## API Access

Due to licensing restrictions, DisGeNET data cannot be bulk downloaded with
a free academic license. Instead, data must be accessed via their API.

### Getting an API Key

1. Apply for academic license: {self.SIGNUP_URL}
2. Approval typically takes 7 business days
3. You'll receive an API key via email

### Setting Up Your API Key

Save your API key in one of these ways:

1. **Environment variable** (recommended):
   ```bash
   export DISGENET_API_KEY=your_key_here
   ```

2. **File** (already done if you provided --disgenet-api-key):
   ```
   echo "your_key_here" > {self.output_dir}/api_key.txt
   ```

## Usage

The DisGeNET client will automatically:
1. Check cache for existing data
2. Query API if not cached (respecting rate limits)
3. Save response to cache for future use

### Cache Structure

```
{cache_dir}/
├── gene_queries/       # Cached gene-disease queries
├── disease_queries/    # Cached disease-gene queries
└── variant_queries/    # Cached variant queries
```

### API Limits

- Free academic license: Reasonable use
- Check https://www.disgenet.org for current limits
- Responses are automatically cached to minimize API calls

## Example Query

```python
from sae_genomics.validation.clients import DisGeNETClient

client = DisGeNETClient(cache_dir="{cache_dir}")
results = client.query_gene("ENSG00000232810")  # TNF gene
```

## Documentation

- API docs: https://www.disgenet.org/api/
- User guide: https://www.disgenet.org/dbinfo

## Version

Cache setup date: {datetime.now().isoformat()}
"""

        with open(self.output_dir / 'README.md', 'w') as f:
            f.write(readme_content)

        print("✓ README.md created with setup instructions")
        print("✓ Cache directory structure created")

        # Create metadata
        metadata = {
            'version': self.get_latest_version(),
            'download_date': datetime.now().isoformat(),
            'database_name': 'disgenet',
            'type': 'api_cache',
            'api_key_configured': bool(self.api_key),
            'cache_dir': str(cache_dir),
            'files': []
        }

        self.save_metadata(metadata)

        if not self.api_key:
            print(f"\n⚠ Next step: Obtain API key from {self.SIGNUP_URL}")

        return metadata

    def verify(self) -> bool:
        """Verify cache structure exists.

        Returns:
            True if cache directories exist
        """
        cache_dir = self.output_dir / 'cache'
        return (
            self.metadata_file.exists() and
            cache_dir.exists() and
            (cache_dir / 'gene_queries').exists()
        )
