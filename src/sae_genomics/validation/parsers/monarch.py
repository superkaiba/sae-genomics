"""Monarch Initiative knowledge graph parser."""

from pathlib import Path
from typing import Iterator, List, Dict, Set, Optional
import duckdb

from .base import GeneDiseaseAssociation, GenePhenotypeAssociation


class MonarchParser:
    """Parser for Monarch Initiative knowledge graph.

    Monarch integrates data from multiple sources (HPO, Orphanet, OMIM, etc.)
    and provides a unified knowledge graph.

    Supports DuckDB format (recommended) or TSV files.
    """

    def __init__(
        self,
        data_dir: Path,
        format: str = 'duckdb',
    ):
        """Initialize Monarch parser.

        Args:
            data_dir: Directory containing Monarch database files
            format: 'duckdb' (recommended) or 'tsv'
        """
        self.data_dir = Path(data_dir)
        self.format = format

        if format == 'duckdb':
            self.db_path = self.data_dir / 'monarch-kg.duckdb'
            if not self.db_path.exists():
                raise FileNotFoundError(f"Monarch DuckDB not found: {self.db_path}")

            # Connect to DuckDB
            self.conn = duckdb.connect(str(self.db_path), read_only=True)

        elif format == 'tsv':
            self.edges_file = self.data_dir / 'monarch-kg-edges.tsv'
            self.nodes_file = self.data_dir / 'monarch-kg-nodes.tsv'

            if not self.edges_file.exists():
                raise FileNotFoundError(f"Monarch edges file not found: {self.edges_file}")

        else:
            raise ValueError(f"Unsupported format: {format}")

    def __del__(self):
        """Close database connection."""
        if hasattr(self, 'conn'):
            self.conn.close()

    def query_gene_diseases(
        self,
        gene_symbol: str,
        include_phenotypes: bool = True,
    ) -> List[GeneDiseaseAssociation]:
        """Query diseases associated with a gene.

        Args:
            gene_symbol: Gene symbol (e.g., "TNF")
            include_phenotypes: Also include phenotype associations

        Returns:
            List of GeneDiseaseAssociation objects
        """
        if self.format != 'duckdb':
            raise NotImplementedError("TSV format not yet supported")

        # Query for gene node (prioritize HGNC/human genes)
        gene_query = """
        SELECT id, name, category
        FROM nodes
        WHERE (name = ? OR id LIKE '%' || ? || '%')
        AND category LIKE '%Gene%'
        ORDER BY CASE WHEN id LIKE 'HGNC:%' THEN 0 ELSE 1 END
        LIMIT 1
        """

        gene_result = self.conn.execute(gene_query, [gene_symbol, gene_symbol]).fetchone()

        if not gene_result:
            return []

        gene_id, gene_name, _ = gene_result

        # Query for disease associations
        disease_query = """
        SELECT DISTINCT
            n2.id as disease_id,
            n2.name as disease_name,
            e.predicate,
            n2.category
        FROM edges e
        JOIN nodes n2 ON e.object = n2.id
        WHERE e.subject = ?
        AND n2.category LIKE '%Disease%'
        """

        results = self.conn.execute(disease_query, [gene_id]).fetchall()

        associations = []
        for disease_id, disease_name, predicate, category in results:
            # Determine evidence level from predicate
            evidence = self._map_predicate_to_evidence(predicate)

            associations.append(GeneDiseaseAssociation(
                gene_id=gene_id,
                gene_symbol=gene_symbol,
                disease_id=disease_id,
                disease_name=disease_name or disease_id,
                source='monarch',
                evidence_level=evidence,
                metadata={
                    'predicate': predicate,
                    'disease_category': category,
                }
            ))

        return associations

    def query_gene_phenotypes(
        self,
        gene_symbol: str,
    ) -> List[GenePhenotypeAssociation]:
        """Query phenotypes associated with a gene.

        Args:
            gene_symbol: Gene symbol

        Returns:
            List of GenePhenotypeAssociation objects
        """
        if self.format != 'duckdb':
            raise NotImplementedError("TSV format not yet supported")

        # Query for gene node (prioritize HGNC/human genes)
        gene_query = """
        SELECT id, name
        FROM nodes
        WHERE (name = ? OR id LIKE '%' || ? || '%')
        AND category LIKE '%Gene%'
        ORDER BY CASE WHEN id LIKE 'HGNC:%' THEN 0 ELSE 1 END
        LIMIT 1
        """

        gene_result = self.conn.execute(gene_query, [gene_symbol, gene_symbol]).fetchone()

        if not gene_result:
            return []

        gene_id, gene_name = gene_result

        # Query for phenotype associations
        phenotype_query = """
        SELECT DISTINCT
            n2.id as phenotype_id,
            n2.name as phenotype_name,
            e.predicate
        FROM edges e
        JOIN nodes n2 ON e.object = n2.id
        WHERE e.subject = ?
        AND n2.category LIKE '%Phenotype%'
        """

        results = self.conn.execute(phenotype_query, [gene_id]).fetchall()

        associations = []
        for phenotype_id, phenotype_name, predicate in results:
            associations.append(GenePhenotypeAssociation(
                gene_id=gene_id,
                gene_symbol=gene_symbol,
                phenotype_id=phenotype_id,
                phenotype_name=phenotype_name or phenotype_id,
                source='monarch',
                metadata={
                    'predicate': predicate,
                }
            ))

        return associations

    def _map_predicate_to_evidence(self, predicate: str) -> str:
        """Map Monarch predicate to evidence level."""
        predicate_lower = predicate.lower()

        if 'causal' in predicate_lower or 'contributes' in predicate_lower:
            return 'definitive'
        elif 'associated' in predicate_lower or 'correlated' in predicate_lower:
            return 'strong'
        elif 'model' in predicate_lower:
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
        for gene in gene_symbols:
            associations = self.query_gene_diseases(gene)
            yield from associations

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
            associations = self.query_gene_diseases(gene)
            diseases = {assoc.disease_id for assoc in associations}
            gene_diseases[gene] = diseases

        return gene_diseases

    def get_all_diseases(self) -> List[Dict[str, str]]:
        """Get all disease nodes from the knowledge graph.

        Returns:
            List of dicts with disease info
        """
        if self.format != 'duckdb':
            raise NotImplementedError("TSV format not yet supported")

        query = """
        SELECT id, name, category
        FROM nodes
        WHERE category LIKE '%Disease%'
        """

        results = self.conn.execute(query).fetchall()

        diseases = []
        for disease_id, name, category in results:
            diseases.append({
                'id': disease_id,
                'name': name,
                'category': category,
            })

        return diseases

    def get_disease_genes(self, disease_id: str) -> List[str]:
        """Get all genes associated with a disease.

        Args:
            disease_id: Disease identifier

        Returns:
            List of gene symbols
        """
        if self.format != 'duckdb':
            raise NotImplementedError("TSV format not yet supported")

        query = """
        SELECT DISTINCT n1.name
        FROM edges e
        JOIN nodes n1 ON e.subject = n1.id
        WHERE e.object = ?
        AND n1.category LIKE '%Gene%'
        """

        results = self.conn.execute(query, [disease_id]).fetchall()

        return [row[0] for row in results if row[0]]
