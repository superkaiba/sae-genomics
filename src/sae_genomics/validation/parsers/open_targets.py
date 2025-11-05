"""Open Targets Platform parser for gene-disease associations."""

from pathlib import Path
from typing import Iterator, List, Dict, Set, Optional
import pandas as pd

from .base import GeneDiseaseAssociation


class OpenTargetsParser:
    """Parser for Open Targets Platform gene-disease associations.

    Open Targets provides comprehensive evidence linking genes to diseases
    from genetics, somatic mutations, drugs, pathways, and literature.

    Data format: Parquet files (columnar format, efficient for queries)
    """

    def __init__(
        self,
        data_dir: Path,
        score_threshold: float = 0.2,
    ):
        """Initialize Open Targets parser.

        Args:
            data_dir: Directory containing Open Targets parquet files
            score_threshold: Minimum association score (0-1)
        """
        self.data_dir = Path(data_dir)
        self.score_threshold = score_threshold

        # Check for required directories
        self.association_dir = self.data_dir / 'association_by_datasource_direct'
        self.disease_dir = self.data_dir / 'disease'
        self.target_dir = self.data_dir / 'target'

        # Load reference data (diseases and targets/genes)
        self._disease_index = None
        self._target_index = None

    def _load_disease_index(self) -> Dict[str, str]:
        """Load disease ID to name mapping."""
        if self._disease_index is not None:
            return self._disease_index

        disease_files = list(self.disease_dir.glob('*.parquet'))

        if not disease_files:
            raise FileNotFoundError(f"No disease parquet files found in {self.disease_dir}")

        # Read all disease files
        dfs = []
        for f in disease_files:
            df = pd.read_parquet(f, columns=['id', 'name'])
            dfs.append(df)

        disease_df = pd.concat(dfs, ignore_index=True)

        # Create index
        self._disease_index = dict(zip(disease_df['id'], disease_df['name']))

        return self._disease_index

    def _load_target_index(self) -> Dict[str, Dict[str, str]]:
        """Load target (gene) ID to symbol/name mapping."""
        if self._target_index is not None:
            return self._target_index

        target_files = list(self.target_dir.glob('*.parquet'))

        if not target_files:
            raise FileNotFoundError(f"No target parquet files found in {self.target_dir}")

        # Read all target files
        dfs = []
        for f in target_files:
            df = pd.read_parquet(f, columns=['id', 'approvedSymbol', 'approvedName'])
            dfs.append(df)

        target_df = pd.concat(dfs, ignore_index=True)

        # Create index
        self._target_index = {}
        for _, row in target_df.iterrows():
            self._target_index[row['id']] = {
                'symbol': row['approvedSymbol'],
                'name': row['approvedName'],
            }

        return self._target_index

    def query_gene(
        self,
        gene_symbol: str,
    ) -> List[GeneDiseaseAssociation]:
        """Query diseases associated with a gene.

        Args:
            gene_symbol: Gene symbol (e.g., "TNF")

        Returns:
            List of GeneDiseaseAssociation objects
        """
        # Load reference data
        target_index = self._load_target_index()
        disease_index = self._load_disease_index()

        # Find target ID for gene symbol
        target_id = None
        for tid, tinfo in target_index.items():
            if tinfo['symbol'] == gene_symbol:
                target_id = tid
                break

        if not target_id:
            return []

        # Load association data
        association_files = list(self.association_dir.glob('*.parquet'))

        if not association_files:
            return []

        associations = []

        # Read association files and filter for this gene
        for assoc_file in association_files:
            # Read parquet file with selected columns
            # Note: Not using filters parameter to ensure compatibility across pandas/pyarrow versions
            df = pd.read_parquet(
                assoc_file,
                columns=['targetId', 'diseaseId', 'score', 'datasourceId']
            )

            # Filter for target gene and score threshold
            df = df[
                (df['targetId'] == target_id) &
                (df['score'] >= self.score_threshold)
            ]

            # Convert to associations
            for _, row in df.iterrows():
                disease_name = disease_index.get(row['diseaseId'], row['diseaseId'])

                associations.append(GeneDiseaseAssociation(
                    gene_id=target_id,
                    gene_symbol=gene_symbol,
                    disease_id=row['diseaseId'],
                    disease_name=disease_name,
                    source=f"open_targets_{row['datasourceId']}",
                    score=float(row['score']),
                    evidence_level=self._map_score_to_evidence(row['score']),
                    metadata={
                        'datasource': row['datasourceId'],
                    }
                ))

        return associations

    def _map_score_to_evidence(self, score: float) -> str:
        """Map association score to evidence level."""
        if score >= 0.7:
            return 'definitive'
        elif score >= 0.5:
            return 'strong'
        elif score >= 0.3:
            return 'moderate'
        else:
            return 'limited'

    def query_genes_batch(
        self,
        gene_symbols: Set[str],
    ) -> Iterator[GeneDiseaseAssociation]:
        """Query multiple genes efficiently.

        Args:
            gene_symbols: Set of gene symbols

        Yields:
            GeneDiseaseAssociation objects
        """
        # Load reference data once
        target_index = self._load_target_index()
        disease_index = self._load_disease_index()

        # Map gene symbols to target IDs
        target_ids = {}
        for tid, tinfo in target_index.items():
            if tinfo['symbol'] in gene_symbols:
                target_ids[tid] = tinfo['symbol']

        if not target_ids:
            return

        # Load all association files
        association_files = list(self.association_dir.glob('*.parquet'))

        for assoc_file in association_files:
            try:
                # Read parquet file (all genes)
                df = pd.read_parquet(
                    assoc_file,
                    columns=['targetId', 'diseaseId', 'score', 'datasourceId']
                )

                # Filter for our genes
                df = df[df['targetId'].isin(target_ids.keys())]
                df = df[df['score'] >= self.score_threshold]

                # Convert to associations
                for _, row in df.iterrows():
                    target_id = row['targetId']
                    gene_symbol = target_ids[target_id]
                    disease_name = disease_index.get(row['diseaseId'], row['diseaseId'])

                    yield GeneDiseaseAssociation(
                        gene_id=target_id,
                        gene_symbol=gene_symbol,
                        disease_id=row['diseaseId'],
                        disease_name=disease_name,
                        source=f"open_targets_{row['datasourceId']}",
                        score=float(row['score']),
                        evidence_level=self._map_score_to_evidence(row['score']),
                        metadata={
                            'datasource': row['datasourceId'],
                        }
                    )

            except Exception:
                continue

    def build_gene_disease_index(
        self,
        gene_symbols: Set[str],
    ) -> Dict[str, Set[str]]:
        """Build index of gene -> diseases.

        Args:
            gene_symbols: Set of gene symbols to index

        Returns:
            Dict mapping gene symbols to sets of disease IDs
        """
        gene_diseases = {}

        for gene in gene_symbols:
            associations = self.query_gene(gene)
            diseases = {assoc.disease_id for assoc in associations}
            gene_diseases[gene] = diseases

        return gene_diseases

    def get_all_diseases(self) -> List[Dict[str, str]]:
        """Get all diseases from Open Targets.

        Returns:
            List of dicts with disease info
        """
        disease_index = self._load_disease_index()

        diseases = []
        for disease_id, name in disease_index.items():
            diseases.append({
                'id': disease_id,
                'name': name,
            })

        return diseases

    def get_disease_genes(
        self,
        disease_id: str,
        score_threshold: Optional[float] = None,
    ) -> List[str]:
        """Get all genes associated with a disease.

        Args:
            disease_id: Disease identifier
            score_threshold: Minimum score (overrides instance threshold)

        Returns:
            List of gene symbols
        """
        threshold = score_threshold or self.score_threshold
        target_index = self._load_target_index()

        genes = set()

        # Load association files
        association_files = list(self.association_dir.glob('*.parquet'))

        for assoc_file in association_files:
            try:
                df = pd.read_parquet(
                    assoc_file,
                    columns=['targetId', 'diseaseId', 'score'],
                    filters=[('diseaseId', '==', disease_id)]
                )

                df = df[df['score'] >= threshold]

                # Map target IDs to gene symbols
                for target_id in df['targetId'].unique():
                    if target_id in target_index:
                        genes.add(target_index[target_id]['symbol'])

            except Exception:
                continue

        return list(genes)
