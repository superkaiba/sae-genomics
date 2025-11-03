"""OMIM API cache setup downloader."""

from pathlib import Path
from typing import List, Dict
from datetime import datetime
import os

from .base import DatabaseDownloader


class OMIMDownloader(DatabaseDownloader):
    """Set up OMIM API cache structure.

    OMIM (Online Mendelian Inheritance in Man) does not allow bulk downloads
    for free academic users. This "downloader" creates the cache structure
    and documentation for API-based access with local caching.

    Database info:
        - Size: API cache (grows with usage)
        - Format: JSON (cached API responses)
        - License: Free for research/education (requires registration)
        - Updates: Real-time via API (daily updates)
        - Limits: 2000 requests/day, API key renewed yearly
        - URL: https://www.omim.org/
    """

    API_BASE = "https://api.omim.org/api"
    SIGNUP_URL = "https://www.omim.org/api"

    def __init__(self, output_dir: Path, api_key: str = None, **kwargs):
        """Initialize OMIM cache setup.

        Args:
            output_dir: Directory for cache structure
            api_key: OMIM API key (optional, can be set later)
            **kwargs: Additional arguments
        """
        super().__init__(output_dir, **kwargs)
        self.api_key = api_key or os.getenv('OMIM_API_KEY')

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
        print(f"\nSetting up OMIM API cache...")
        print(f"Output directory: {self.output_dir}")

        # Create cache directories
        cache_dir = self.output_dir / 'cache'
        cache_dir.mkdir(parents=True, exist_ok=True)

        (cache_dir / 'gene_pheno').mkdir(exist_ok=True)
        (cache_dir / 'morbid_map').mkdir(exist_ok=True)
        (cache_dir / 'clinical_synopsis').mkdir(exist_ok=True)

        # Save API key if provided
        if self.api_key:
            with open(self.output_dir / 'api_key.txt', 'w') as f:
                f.write(self.api_key)
            print("✓ API key saved")
        else:
            print("⚠ No API key provided. You'll need to obtain one from:")
            print(f"  {self.SIGNUP_URL}")

        # Create README with instructions
        readme_content = f"""# OMIM API Cache

## About OMIM

OMIM (Online Mendelian Inheritance in Man) is a comprehensive database of
human genes and genetic phenotypes, focusing on Mendelian disorders.

## API Access

OMIM requires API access for all programmatic queries. Bulk downloads are
restricted to paid commercial licenses.

### Getting an API Key

1. Register at: {self.SIGNUP_URL}
2. Fill out the API key request form
3. **Approval is typically within 2 hours**
4. You'll receive an API key via email
5. **Key must be renewed annually**

### API Limits

- **2000 requests per day** (resets at midnight UTC)
- Responses are cached locally to minimize API usage
- Plan your queries accordingly

### Setting Up Your API Key

Save your API key in one of these ways:

1. **Environment variable** (recommended):
   ```bash
   export OMIM_API_KEY=your_key_here
   ```

2. **File** (already done if you provided --omim-api-key):
   ```
   echo "your_key_here" > {self.output_dir}/api_key.txt
   ```

## Usage

The OMIM client will automatically:
1. Check cache for existing data (30-day TTL)
2. Query API if not cached (respecting daily limits)
3. Save response to cache for future use

### Cache Structure

```
{cache_dir}/
├── gene_pheno/         # Cached gene-phenotype queries
├── morbid_map/         # Cached morbid map queries
└── clinical_synopsis/  # Cached clinical synopsis data
```

### Cache Policy

- **TTL**: 30 days (longer than other databases due to stability)
- **Format**: JSON files named by MIM number or gene ID
- **Automatic cleanup**: Old cache files can be removed after TTL

## Example Query

```python
from sae_genomics.validation.clients import OMIMClient

client = OMIMClient(cache_dir="{cache_dir}")
results = client.query_gene_symbol("TNF")
# Returns: Gene-phenotype associations for TNF
```

## Response Format

OMIM API returns JSON with:
- MIM numbers (6-digit identifiers)
- Gene-phenotype relationships
- Clinical descriptions
- Allelic variants
- References

## Documentation

- API docs: https://www.omim.org/api
- Data structure: https://www.omim.org/help/api
- Search guide: https://www.omim.org/help/search

## Important Notes

1. **Daily limit**: Plan batch queries carefully
2. **Cache aggressively**: OMIM data is relatively stable
3. **Renewal**: Remember to renew API key annually
4. **Attribution**: Cite OMIM in any publications

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
            'database_name': 'omim',
            'type': 'api_cache',
            'api_key_configured': bool(self.api_key),
            'cache_dir': str(cache_dir),
            'daily_limit': 2000,
            'cache_ttl_days': 30,
            'files': []
        }

        self.save_metadata(metadata)

        if not self.api_key:
            print(f"\n⚠ Next step: Obtain API key from {self.SIGNUP_URL}")
            print("  (Approval typically within 2 hours)")

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
            (cache_dir / 'gene_pheno').exists()
        )
