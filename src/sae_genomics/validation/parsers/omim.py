"""OMIM API client and parser for gene-disease associations."""

from pathlib import Path
from typing import Iterator, List, Dict, Optional, Set
import json
import time
import requests
import hashlib

from .base import GeneDiseaseAssociation


class OMIMParser:
    """API client for OMIM (Online Mendelian Inheritance in Man).

    OMIM requires an API key for access. Free for research/educational use.
    API documentation: https://www.omim.org/api

    Rate limits: ~2000 requests per day for free API keys.
    """

    API_BASE = "https://api.omim.org/api"

    def __init__(
        self,
        data_dir: Path,
        api_key: Optional[str] = None,
        cache_ttl_days: int = 30,
        daily_limit: int = 2000,
    ):
        """Initialize OMIM parser.

        Args:
            data_dir: Directory containing cache structure
            api_key: OMIM API key (from environment or file)
            cache_ttl_days: Cache validity in days
            daily_limit: Daily API request limit
        """
        self.data_dir = Path(data_dir)
        self.cache_ttl_days = cache_ttl_days
        self.daily_limit = daily_limit

        # Load API key
        self.api_key = api_key or self._load_api_key()

        # Cache directories
        self.cache_dir = self.data_dir / 'cache'
        self.gene_cache = self.cache_dir / 'gene_queries'
        self.entry_cache = self.cache_dir / 'entry_queries'

        # Create cache dirs
        self.gene_cache.mkdir(parents=True, exist_ok=True)
        self.entry_cache.mkdir(parents=True, exist_ok=True)

        # Track API usage
        self.usage_file = self.cache_dir / 'api_usage.json'

    def _load_api_key(self) -> Optional[str]:
        """Load API key from file or environment."""
        import os

        # Try environment variable
        api_key = os.getenv('OMIM_API_KEY')
        if api_key:
            return api_key

        # Try api_key.txt file
        key_file = self.data_dir / 'api_key.txt'
        if key_file.exists():
            return key_file.read_text().strip()

        return None

    def _get_cache_path(self, query_type: str, identifier: str) -> Path:
        """Get cache file path for a query."""
        id_hash = hashlib.md5(identifier.encode()).hexdigest()

        if query_type == 'gene':
            return self.gene_cache / f"{id_hash}.json"
        elif query_type == 'entry':
            return self.entry_cache / f"{id_hash}.json"
        else:
            raise ValueError(f"Unknown query type: {query_type}")

    def _is_cache_valid(self, cache_path: Path) -> bool:
        """Check if cache file is still valid."""
        if not cache_path.exists():
            return False

        age_days = (time.time() - cache_path.stat().st_mtime) / 86400
        return age_days < self.cache_ttl_days

    def _check_rate_limit(self) -> bool:
        """Check if we're within daily rate limit."""
        if not self.usage_file.exists():
            return True

        with open(self.usage_file) as f:
            usage = json.load(f)

        today = time.strftime('%Y-%m-%d')
        daily_count = usage.get(today, 0)

        return daily_count < self.daily_limit

    def _increment_usage(self):
        """Increment API usage counter."""
        if self.usage_file.exists():
            with open(self.usage_file) as f:
                usage = json.load(f)
        else:
            usage = {}

        today = time.strftime('%Y-%m-%d')
        usage[today] = usage.get(today, 0) + 1

        with open(self.usage_file, 'w') as f:
            json.dump(usage, f)

    def _query_api(self, endpoint: str, params: Dict) -> Dict:
        """Query OMIM API with rate limiting."""
        if not self.api_key:
            raise ValueError(
                "OMIM API key required. Get one at https://www.omim.org/api\n"
                "Set OMIM_API_KEY environment variable or save to api_key.txt"
            )

        if not self._check_rate_limit():
            raise ValueError(f"Daily API limit ({self.daily_limit}) reached")

        # Add API key to params
        params['apiKey'] = self.api_key
        params['format'] = 'json'

        url = f"{self.API_BASE}{endpoint}"

        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()

            self._increment_usage()
            time.sleep(0.5)  # Rate limiting

            return response.json()
        except requests.HTTPError as e:
            if e.response.status_code == 401:
                raise ValueError("Invalid API key") from e
            elif e.response.status_code == 429:
                raise ValueError("API rate limit exceeded") from e
            else:
                raise

    def query_gene(
        self,
        gene_symbol: str,
        use_cache: bool = True
    ) -> List[GeneDiseaseAssociation]:
        """Query gene-disease associations for a gene.

        Args:
            gene_symbol: Gene symbol (e.g., "TNF")
            use_cache: Whether to use cached results

        Returns:
            List of GeneDiseaseAssociation objects
        """
        cache_path = self._get_cache_path('gene', gene_symbol)

        # Check cache
        if use_cache and self._is_cache_valid(cache_path):
            with open(cache_path) as f:
                data = json.load(f)
        else:
            # Query API
            if not self.api_key:
                return []

            try:
                # Search for gene
                data = self._query_api('/entry/search', {
                    'search': f'gene_symbol:{gene_symbol}',
                    'include': 'geneMap',
                    'retrieve': 'geneMap',
                })

                # Cache response
                with open(cache_path, 'w') as f:
                    json.dump(data, f)

            except Exception:
                return []

        # Parse results
        associations = []

        search_response = data.get('omim', {}).get('searchResponse', {})
        entry_list = search_response.get('entryList', [])

        for entry in entry_list:
            entry_data = entry.get('entry', {})

            # Get gene map info
            gene_map = entry_data.get('geneMap', {})
            phenotype_map_list = gene_map.get('phenotypeMapList', [])

            for pheno_map in phenotype_map_list:
                pheno = pheno_map.get('phenotypeMap', {})

                # Get phenotype MIM number and name
                mim_number = pheno.get('mimNumber')
                phenotype = pheno.get('phenotype', '')

                if not mim_number or not phenotype:
                    continue

                # Parse phenotype name and inheritance
                # Format is usually: "Phenotype name, inheritance"
                parts = phenotype.split(',')
                disease_name = parts[0].strip() if parts else phenotype

                # Get mapping key (indicates association type)
                mapping_key = pheno.get('phenotypeMappingKey', '')

                associations.append(GeneDiseaseAssociation(
                    gene_id=str(entry_data.get('mimNumber', '')),
                    gene_symbol=gene_symbol,
                    disease_id=f"OMIM:{mim_number}",
                    disease_name=disease_name,
                    source='omim',
                    evidence_level=self._map_phenotype_key(mapping_key),
                    metadata={
                        'gene_mim': entry_data.get('mimNumber'),
                        'phenotype_mim': mim_number,
                        'mapping_key': mapping_key,
                        'inheritance': ','.join(parts[1:]) if len(parts) > 1 else '',
                        'gene_symbols': gene_map.get('geneSymbols', ''),
                    }
                ))

        return associations

    def _map_phenotype_key(self, key: str) -> str:
        """Map OMIM phenotype mapping key to evidence level.

        Mapping keys:
        1 - disorder was positioned by mapping of the wildtype gene
        2 - disorder was positioned by deletion or duplication mapping
        3 - molecular basis of the disorder is known
        4 - disorder is a chromosome deletion or duplication syndrome
        """
        key_map = {
            '1': 'limited',
            '2': 'moderate',
            '3': 'definitive',
            '4': 'limited',
        }
        return key_map.get(str(key), 'unknown')

    def query_genes_batch(
        self,
        gene_symbols: Set[str],
        use_cache: bool = True
    ) -> Iterator[GeneDiseaseAssociation]:
        """Query multiple genes efficiently.

        Args:
            gene_symbols: Set of gene symbols
            use_cache: Whether to use cached results

        Yields:
            GeneDiseaseAssociation objects
        """
        for gene in gene_symbols:
            associations = self.query_gene(gene, use_cache=use_cache)
            yield from associations

    def build_gene_disease_index(
        self,
        gene_symbols: Set[str],
        use_cache: bool = True
    ) -> Dict[str, Set[str]]:
        """Build index of gene -> diseases.

        Args:
            gene_symbols: Set of gene symbols to index
            use_cache: Whether to use cached results

        Returns:
            Dict mapping gene symbols to sets of disease IDs
        """
        gene_diseases = {}

        for gene in gene_symbols:
            associations = self.query_gene(gene, use_cache=use_cache)
            diseases = {assoc.disease_id for assoc in associations}
            gene_diseases[gene] = diseases

        return gene_diseases
