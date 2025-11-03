"""Base classes for database parsers."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterator, List, Optional, Dict, Set
from dataclasses import dataclass, field, asdict
import json


@dataclass
class GeneDiseaseAssociation:
    """Standardized gene-disease association.

    This represents a link between a gene and a disease, with metadata
    about the evidence supporting the association.
    """

    gene_id: str  # Ensembl gene ID (standardized format: ENSG...)
    gene_symbol: str  # HGNC gene symbol (e.g., "TNF", "IL6")
    disease_id: str  # Disease ID (preferably MONDO, EFO, or OMIM)
    disease_name: str  # Human-readable disease name
    source: str  # Database name (e.g., "gencc", "orphanet")

    # Optional fields
    evidence_level: Optional[str] = None  # e.g., "definitive", "strong", "moderate"
    score: Optional[float] = None  # Confidence score (0-1) if available
    pmids: List[str] = field(default_factory=list)  # PubMed IDs supporting association
    metadata: Dict = field(default_factory=dict)  # Additional database-specific info

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return asdict(self)

    def to_json(self) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict())


@dataclass
class GenePhenotypeAssociation:
    """Standardized gene-phenotype association.

    Similar to GeneDiseaseAssociation but for phenotypes (HPO terms).
    """

    gene_id: str  # Ensembl gene ID
    gene_symbol: str  # HGNC gene symbol
    phenotype_id: str  # HPO ID (e.g., "HP:0001234")
    phenotype_name: str  # Human-readable phenotype name
    source: str  # Database name

    # Optional fields
    evidence_code: Optional[str] = None  # e.g., "IEA", "TAS", "PCS"
    frequency: Optional[str] = None  # e.g., "Frequent", "Occasional"
    onset: Optional[str] = None  # Age of onset
    pmids: List[str] = field(default_factory=list)
    metadata: Dict = field(default_factory=dict)

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return asdict(self)


@dataclass
class GeneGeneInteraction:
    """Protein-protein interaction from STRING."""

    gene_id_1: str  # First gene (Ensembl ID)
    gene_id_2: str  # Second gene (Ensembl ID)
    combined_score: float  # Combined confidence score (0-1000)
    source: str = "string"

    # Individual scores by evidence type
    experimental_score: Optional[float] = None
    database_score: Optional[float] = None
    text_mining_score: Optional[float] = None
    coexpression_score: Optional[float] = None

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return asdict(self)


class DatabaseParser(ABC):
    """Base class for database parsers.

    Parsers read downloaded database files and yield standardized
    association objects.
    """

    def __init__(self, data_dir: Path):
        """Initialize parser.

        Args:
            data_dir: Directory containing downloaded database files
        """
        self.data_dir = Path(data_dir)

        if not self.data_dir.exists():
            raise FileNotFoundError(f"Database directory not found: {self.data_dir}")

    @abstractmethod
    def parse_gene_disease_associations(self) -> Iterator[GeneDiseaseAssociation]:
        """Parse and yield gene-disease associations.

        Yields:
            GeneDiseaseAssociation objects

        Raises:
            NotImplementedError: If the database doesn't provide disease associations
        """
        pass

    def parse_gene_phenotype_associations(self) -> Iterator[GenePhenotypeAssociation]:
        """Parse and yield gene-phenotype associations.

        Yields:
            GenePhenotypeAssociation objects

        Raises:
            NotImplementedError: If the database doesn't provide phenotype associations
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} does not provide phenotype associations"
        )

    def parse_gene_interactions(self) -> Iterator[GeneGeneInteraction]:
        """Parse and yield gene-gene interactions.

        Yields:
            GeneGeneInteraction objects

        Raises:
            NotImplementedError: If the database doesn't provide interactions
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} does not provide gene interactions"
        )

    def get_all_genes(self) -> Set[str]:
        """Get all unique gene IDs in the database.

        Returns:
            Set of gene IDs (Ensembl format)
        """
        genes = set()
        try:
            for assoc in self.parse_gene_disease_associations():
                genes.add(assoc.gene_id)
        except NotImplementedError:
            pass

        try:
            for assoc in self.parse_gene_phenotype_associations():
                genes.add(assoc.gene_id)
        except NotImplementedError:
            pass

        return genes

    def get_all_diseases(self) -> Set[str]:
        """Get all unique disease IDs in the database.

        Returns:
            Set of disease IDs
        """
        diseases = set()
        try:
            for assoc in self.parse_gene_disease_associations():
                diseases.add(assoc.disease_id)
        except NotImplementedError:
            pass

        return diseases

    def count_associations(self) -> Dict[str, int]:
        """Count associations by type.

        Returns:
            Dict with counts for each association type
        """
        counts = {}

        try:
            gene_disease = sum(1 for _ in self.parse_gene_disease_associations())
            counts['gene_disease'] = gene_disease
        except NotImplementedError:
            pass

        try:
            gene_phenotype = sum(1 for _ in self.parse_gene_phenotype_associations())
            counts['gene_phenotype'] = gene_phenotype
        except NotImplementedError:
            pass

        try:
            gene_interactions = sum(1 for _ in self.parse_gene_interactions())
            counts['gene_interactions'] = gene_interactions
        except NotImplementedError:
            pass

        return counts
