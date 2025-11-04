"""DisGeNET API client and parser for gene-disease associations."""

from pathlib import Path
from typing import Iterator, List, Dict, Optional, Set
import json
import time
import requests
import hashlib

from .base import GeneDiseaseAssociation


class DisGeNETParser:
    """API client for DisGeNET gene-disease associations.

    DisGeNET doesn't provide bulk downloads for free academic users.
    This parser queries the API and caches responses locally.

    API documentation: https://www.disgenet.org/api/
    """

    API_BASE = "https://www.disgenet.org/api"

    def __init__(
        self,
        data_dir: Path,
        api_key: Optional[str] = None,
        score_threshold: float = 0.3,
        cache_ttl_days: int = 30,
    ):
        """Initialize DisGeNET parser.

        Args:
            data_dir: Directory containing cache structure
            api_key: DisGeNET API key (from environment or file)
            score_threshold: Minimum association score (0-1)
            cache_ttl_days: Cache validity in days
        """
        self.data_dir = Path(data_dir)
        self.score_threshold = score_threshold
        self.cache_ttl_days = cache_ttl_days

        # Load API key
        self.api_key = api_key or self._load_api_key()

        # Cache directories
        self.cache_dir = self.data_dir / 'cache'
        self.gene_cache = self.cache_dir / 'gene_queries'
        self.disease_cache = self.cache_dir / 'disease_queries'

        # Create cache dirs if needed
        self.gene_cache.mkdir(parents=True, exist_ok=True)
        self.disease_cache.mkdir(parents=True, exist_ok=True)

    def _load_api_key(self) -> Optional[str]:
        """Load API key from file or environment."""
        import os

        # Try environment variable first
        api_key = os.getenv('DISGENET_API_KEY')
        if api_key:
            return api_key

        # Try api_key.txt file
        key_file = self.data_dir / 'api_key.txt'
        if key_file.exists():
            return key_file.read_text().strip()

        return None

    def _get_cache_path(self, query_type: str, identifier: str) -> Path:
        """Get cache file path for a query."""
        # Use hash to avoid filesystem issues
        id_hash = hashlib.md5(identifier.encode()).hexdigest()

        if query_type == 'gene':
            return self.gene_cache / f"{id_hash}.json"
        elif query_type == 'disease':
            return self.disease_cache / f"{id_hash}.json"
        else:
            raise ValueError(f"Unknown query type: {query_type}")

    def _is_cache_valid(self, cache_path: Path) -> bool:
        """Check if cache file is still valid."""
        if not cache_path.exists():
            return False

        # Check age
        age_days = (time.time() - cache_path.stat().st_mtime) / 86400
        return age_days < self.cache_ttl_days

    def _query_api(self, endpoint: str, params: Dict) -> Dict:
        """Query DisGeNET API with rate limiting."""
        if not self.api_key:
            raise ValueError(
                "DisGeNET API key required. Set DISGENET_API_KEY environment "
                "variable or save key to api_key.txt"
            )

        headers = {
            'Authorization': f'Bearer {self.api_key}'
        }

        url = f"{self.API_BASE}{endpoint}"

        try:
            response = requests.get(url, headers=headers, params=params, timeout=30)
            response.raise_for_status()
            return response.json()
        except requests.HTTPError as e:
            if e.response.status_code == 401:
                raise ValueError("Invalid API key") from e
            elif e.response.status_code == 429:
                raise ValueError("API rate limit exceeded") from e
            else:
                raise

    def query_gene(self, gene_symbol: str, use_cache: bool = True) -> List[GeneDiseaseAssociation]:
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
                return []  # Return empty if no API key

            try:
                data = self._query_api('/gda/gene', {'gene': gene_symbol})

                # Cache response
                with open(cache_path, 'w') as f:
                    json.dump(data, f)

                # Rate limiting
                time.sleep(0.5)
            except Exception:
                # Return empty on API errors
                return []

        # Parse results
        associations = []
        for item in data.get('results', []):
            score = float(item.get('score', 0))
            if score < self.score_threshold:
                continue

            associations.append(GeneDiseaseAssociation(
                gene_id=item.get('geneid', ''),
                gene_symbol=gene_symbol,
                disease_id=item.get('diseaseid', ''),
                disease_name=item.get('disease_name', ''),
                source='disgenet',
                score=score,
                evidence_level=item.get('el', ''),
                pmids=item.get('pmids', '').split(';') if item.get('pmids') else [],
                metadata={
                    'disease_type': item.get('disease_type'),
                    'disease_class': item.get('disease_class'),
                    'n_pmids': item.get('NofPmids', 0),
                    'n_snps': item.get('NofSnps', 0),
                }
            ))

        return associations

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

    def get_disease_genes(
        self,
        disease_id: str,
        use_cache: bool = True
    ) -> List[GeneDiseaseAssociation]:
        """Get all genes associated with a disease.

        Args:
            disease_id: Disease ID (e.g., "C0006142" for breast cancer)
            use_cache: Whether to use cached results

        Returns:
            List of GeneDiseaseAssociation objects
        """
        cache_path = self._get_cache_path('disease', disease_id)

        # Check cache
        if use_cache and self._is_cache_valid(cache_path):
            with open(cache_path) as f:
                data = json.load(f)
        else:
            # Query API
            if not self.api_key:
                return []

            try:
                data = self._query_api('/gda/disease', {'disease': disease_id})

                # Cache response
                with open(cache_path, 'w') as f:
                    json.dump(data, f)

                # Rate limiting
                time.sleep(0.5)
            except Exception:
                return []

        # Parse results
        associations = []
        for item in data.get('results', []):
            score = float(item.get('score', 0))
            if score < self.score_threshold:
                continue

            associations.append(GeneDiseaseAssociation(
                gene_id=item.get('geneid', ''),
                gene_symbol=item.get('gene_symbol', ''),
                disease_id=disease_id,
                disease_name=item.get('disease_name', ''),
                source='disgenet',
                score=score,
                evidence_level=item.get('el', ''),
                pmids=item.get('pmids', '').split(';') if item.get('pmids') else [],
                metadata={
                    'disease_type': item.get('disease_type'),
                    'n_pmids': item.get('NofPmids', 0),
                    'n_snps': item.get('NofSnps', 0),
                }
            ))

        return associations

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
