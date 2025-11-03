"""STRING database parser for protein-protein interactions."""

from pathlib import Path
from typing import Iterator, Dict, Set, Tuple, Optional
import csv


class ProteinInteraction:
    """Represents a protein-protein interaction from STRING."""

    def __init__(
        self,
        protein_a: str,
        protein_b: str,
        combined_score: float,
        experimental_score: float = 0.0,
        database_score: float = 0.0,
        textmining_score: float = 0.0,
        coexpression_score: float = 0.0,
    ):
        """Initialize protein interaction.

        Args:
            protein_a: First protein (STRING ID or gene symbol)
            protein_b: Second protein (STRING ID or gene symbol)
            combined_score: Combined confidence score (0-1000)
            experimental_score: Experimental evidence score
            database_score: Database evidence score
            textmining_score: Text mining evidence score
            coexpression_score: Co-expression evidence score
        """
        self.protein_a = protein_a
        self.protein_b = protein_b
        self.combined_score = combined_score
        self.experimental_score = experimental_score
        self.database_score = database_score
        self.textmining_score = textmining_score
        self.coexpression_score = coexpression_score

    @property
    def confidence(self) -> float:
        """Get confidence score normalized to 0-1."""
        return self.combined_score / 1000.0

    def __repr__(self) -> str:
        return f"ProteinInteraction({self.protein_a} - {self.protein_b}, score={self.combined_score})"


class STRINGParser:
    """Parser for STRING protein-protein interaction network.

    STRING provides known and predicted protein-protein interactions
    with confidence scores based on multiple evidence types.

    File formats:
        - protein.links.detailed.txt: Detailed interaction scores
        - protein.links.txt: Simple interaction scores
        - protein.info.txt: Protein metadata (mapping to gene symbols)
    """

    def __init__(
        self,
        data_dir: Path,
        score_threshold: int = 400,
        organism: str = '9606',  # Human
    ):
        """Initialize STRING parser.

        Args:
            data_dir: Directory containing STRING files
            score_threshold: Minimum combined score (0-1000, default: 400 = medium confidence)
            organism: NCBI taxonomy ID (default: 9606 for human)
        """
        self.data_dir = data_dir
        self.organism = organism
        self.score_threshold = score_threshold

        # File paths
        self.links_detailed_file = self.data_dir / f'{organism}.protein.links.detailed.v12.0.txt'
        self.links_file = self.data_dir / f'{organism}.protein.links.v12.0.txt'
        self.info_file = self.data_dir / f'{organism}.protein.info.v12.0.txt'

        # Check which files exist
        if not self.links_detailed_file.exists() and not self.links_file.exists():
            raise FileNotFoundError(
                f"STRING files not found in {self.data_dir}. "
                f"Expected: {self.links_detailed_file.name} or {self.links_file.name}"
            )

        # Load protein ID to gene symbol mapping
        self._protein_to_gene: Optional[Dict[str, str]] = None

    def _load_protein_info(self) -> Dict[str, str]:
        """Load mapping from STRING protein IDs to gene symbols.

        Returns:
            Dict mapping protein IDs to gene symbols
        """
        if self._protein_to_gene is not None:
            return self._protein_to_gene

        self._protein_to_gene = {}

        if not self.info_file.exists():
            return self._protein_to_gene

        with open(self.info_file, 'r', encoding='utf-8') as f:
            # Skip header
            header = f.readline()

            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue

                protein_id = parts[0]  # e.g., "9606.ENSP00000000233"
                preferred_name = parts[1]  # Gene symbol

                self._protein_to_gene[protein_id] = preferred_name

        return self._protein_to_gene

    def _protein_id_to_gene(self, protein_id: str) -> str:
        """Convert STRING protein ID to gene symbol.

        Args:
            protein_id: STRING protein ID (e.g., "9606.ENSP00000000233")

        Returns:
            Gene symbol or protein ID if mapping not available
        """
        mapping = self._load_protein_info()
        return mapping.get(protein_id, protein_id)

    def parse_interactions(self, detailed: bool = True) -> Iterator[ProteinInteraction]:
        """Parse STRING protein-protein interactions.

        Args:
            detailed: Whether to parse detailed scores (default: True)

        Yields:
            ProteinInteraction objects
        """
        if detailed and self.links_detailed_file.exists():
            yield from self._parse_detailed_links()
        elif self.links_file.exists():
            yield from self._parse_simple_links()
        else:
            raise FileNotFoundError("No STRING interaction files found")

    def _parse_detailed_links(self) -> Iterator[ProteinInteraction]:
        """Parse detailed STRING interactions with evidence scores."""
        with open(self.links_detailed_file, 'r', encoding='utf-8') as f:
            # Skip header
            header = f.readline()

            for line in f:
                parts = line.strip().split()
                if len(parts) < 10:
                    continue

                protein1 = parts[0]
                protein2 = parts[1]
                combined_score = int(parts[9])

                # Filter by score threshold
                if combined_score < self.score_threshold:
                    continue

                # Convert to gene symbols
                gene1 = self._protein_id_to_gene(protein1)
                gene2 = self._protein_id_to_gene(protein2)

                yield ProteinInteraction(
                    protein_a=gene1,
                    protein_b=gene2,
                    combined_score=combined_score,
                    experimental_score=int(parts[5]) if len(parts) > 5 else 0,
                    database_score=int(parts[6]) if len(parts) > 6 else 0,
                    textmining_score=int(parts[7]) if len(parts) > 7 else 0,
                    coexpression_score=int(parts[4]) if len(parts) > 4 else 0,
                )

    def _parse_simple_links(self) -> Iterator[ProteinInteraction]:
        """Parse simple STRING interactions (combined score only)."""
        with open(self.links_file, 'r', encoding='utf-8') as f:
            # Skip header
            header = f.readline()

            for line in f:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue

                protein1 = parts[0]
                protein2 = parts[1]
                combined_score = int(parts[2])

                # Filter by score threshold
                if combined_score < self.score_threshold:
                    continue

                # Convert to gene symbols
                gene1 = self._protein_id_to_gene(protein1)
                gene2 = self._protein_id_to_gene(protein2)

                yield ProteinInteraction(
                    protein_a=gene1,
                    protein_b=gene2,
                    combined_score=combined_score,
                )

    def get_interactors(
        self,
        gene_symbol: str,
        min_score: Optional[int] = None,
    ) -> Set[str]:
        """Get all proteins that interact with a given gene.

        Args:
            gene_symbol: Gene symbol to query
            min_score: Minimum interaction score (default: use instance threshold)

        Returns:
            Set of interacting gene symbols
        """
        min_score = min_score or self.score_threshold
        interactors = set()

        for interaction in self.parse_interactions():
            if interaction.combined_score < min_score:
                continue

            if interaction.protein_a.upper() == gene_symbol.upper():
                interactors.add(interaction.protein_b)
            elif interaction.protein_b.upper() == gene_symbol.upper():
                interactors.add(interaction.protein_a)

        return interactors

    def get_interaction_score(self, gene_a: str, gene_b: str) -> Optional[float]:
        """Get interaction score between two genes.

        Args:
            gene_a: First gene symbol
            gene_b: Second gene symbol

        Returns:
            Combined score (0-1000) or None if no interaction
        """
        for interaction in self.parse_interactions():
            if ((interaction.protein_a.upper() == gene_a.upper() and
                 interaction.protein_b.upper() == gene_b.upper()) or
                (interaction.protein_a.upper() == gene_b.upper() and
                 interaction.protein_b.upper() == gene_a.upper())):
                return interaction.combined_score

        return None

    def build_network(
        self,
        gene_set: Set[str],
        min_score: Optional[int] = None,
    ) -> Dict[str, Set[str]]:
        """Build interaction network for a set of genes.

        Args:
            gene_set: Set of gene symbols
            min_score: Minimum interaction score

        Returns:
            Dict mapping each gene to its interactors within the gene set
        """
        min_score = min_score or self.score_threshold
        network = {gene: set() for gene in gene_set}
        gene_set_upper = {g.upper() for g in gene_set}

        for interaction in self.parse_interactions():
            if interaction.combined_score < min_score:
                continue

            gene_a = interaction.protein_a
            gene_b = interaction.protein_b

            # Check if both genes are in the query set
            if gene_a.upper() in gene_set_upper and gene_b.upper() in gene_set_upper:
                # Find original case in gene_set
                gene_a_original = next((g for g in gene_set if g.upper() == gene_a.upper()), gene_a)
                gene_b_original = next((g for g in gene_set if g.upper() == gene_b.upper()), gene_b)

                network[gene_a_original].add(gene_b_original)
                network[gene_b_original].add(gene_a_original)

        return network
