"""HPO (Human Phenotype Ontology) database parser."""

from pathlib import Path
from typing import Iterator, Dict, Optional
import csv

from .base import DatabaseParser, GenePhenotypeAssociation, GeneDiseaseAssociation


class HPOParser(DatabaseParser):
    """Parser for HPO gene-phenotype associations.

    HPO provides:
    1. Ontology structure (hp-base.obo)
    2. Gene-phenotype associations (genes_to_phenotype.txt)
    3. Gene-disease associations via phenotypes (genes_to_disease.txt)

    File formats: TSV files
    """

    def __init__(self, data_dir: Path, load_ontology: bool = False):
        """Initialize HPO parser.

        Args:
            data_dir: Directory containing HPO files
            load_ontology: Whether to load the OBO ontology file (requires pronto)
        """
        super().__init__(data_dir)

        self.genes_to_phenotype_file = self.data_dir / 'genes_to_phenotype.txt'
        self.genes_to_disease_file = self.data_dir / 'genes_to_disease.txt'
        self.obo_file = self.data_dir / 'hp-base.obo'

        # Check required files exist
        if not self.genes_to_phenotype_file.exists():
            raise FileNotFoundError(
                f"HPO genes_to_phenotype file not found: {self.genes_to_phenotype_file}"
            )

        self.ontology = None
        if load_ontology:
            self._load_ontology()

    def _load_ontology(self):
        """Load HPO ontology using pronto."""
        try:
            import pronto
            if self.obo_file.exists():
                self.ontology = pronto.Ontology(str(self.obo_file))
            else:
                print(f"Warning: OBO file not found: {self.obo_file}")
        except ImportError:
            print("Warning: pronto library not available. Ontology features disabled.")

    def parse_gene_phenotype_associations(self) -> Iterator[GenePhenotypeAssociation]:
        """Parse HPO gene-phenotype associations.

        File format (genes_to_phenotype.txt):
            ncbi_gene_id    gene_symbol    hpo_id    hpo_name    frequency    ...

        Yields:
            GenePhenotypeAssociation objects
        """
        with open(self.genes_to_phenotype_file, 'r', encoding='utf-8') as f:
            # Process all lines - check each line to see if it's a header
            for line in f:
                # Skip comment lines and header lines
                if line.startswith('#') or 'ncbi_gene_id' in line.lower() or 'gene_symbol' in line.lower():
                    continue

                parts = line.strip().split('\t')

                if len(parts) < 4:
                    continue

                ncbi_gene_id = parts[0]
                gene_symbol = parts[1]
                hpo_id = parts[2]
                hpo_name = parts[3]
                frequency = parts[4] if len(parts) > 4 else None
                disease_id = parts[5] if len(parts) > 5 else None
                qualifier = parts[6] if len(parts) > 6 else None
                aspect = parts[7] if len(parts) > 7 else None
                evidence = parts[9] if len(parts) > 9 else None

                if qualifier and qualifier.strip() == 'NOT':
                    # Skip negative associations
                    continue

                yield GenePhenotypeAssociation(
                    gene_id=ncbi_gene_id,
                    gene_symbol=gene_symbol,
                    phenotype_id=hpo_id,
                    phenotype_name=hpo_name,
                    source='hpo',
                    frequency=frequency,
                    evidence_code=evidence,
                    metadata={
                        'disease_id': disease_id,
                        'aspect': aspect,
                        'qualifier': qualifier,
                    }
                )

    def parse_gene_disease_associations(self) -> Iterator[GeneDiseaseAssociation]:
        """Parse HPO gene-disease associations.

        File format (genes_to_disease.txt):
            ncbi_gene_id    gene_symbol    disease_id    disease_name    ...

        Yields:
            GeneDiseaseAssociation objects
        """
        if not self.genes_to_disease_file.exists():
            return

        with open(self.genes_to_disease_file, 'r', encoding='utf-8') as f:
            # Process all lines - check each line to see if it's a header
            for line in f:
                # Skip comment lines and header lines
                if line.startswith('#') or 'ncbi_gene_id' in line.lower() or 'gene_symbol' in line.lower():
                    continue
                parts = line.strip().split('\t')

                if len(parts) < 4:
                    continue

                ncbi_gene_id = parts[0]
                gene_symbol = parts[1]
                disease_id = parts[2]  # OMIM ID or ORPHA ID
                disease_name = parts[3]

                yield GeneDiseaseAssociation(
                    gene_id=ncbi_gene_id,  # Will be mapped to Ensembl later
                    gene_symbol=gene_symbol,
                    disease_id=disease_id,
                    disease_name=disease_name,
                    source='hpo',
                    evidence_level='curated',  # HPO is manually curated
                )

    def get_term_ancestors(self, hpo_id: str) -> list:
        """Get ancestor terms for an HPO term.

        Requires ontology to be loaded.

        Args:
            hpo_id: HPO ID (e.g., "HP:0001234")

        Returns:
            List of ancestor HPO IDs
        """
        if not self.ontology:
            return []

        try:
            term = self.ontology[hpo_id]
            return [str(anc.id) for anc in term.superclasses() if str(anc.id) != hpo_id]
        except KeyError:
            return []

    def get_term_descendants(self, hpo_id: str) -> list:
        """Get descendant terms for an HPO term.

        Args:
            hpo_id: HPO ID

        Returns:
            List of descendant HPO IDs
        """
        if not self.ontology:
            return []

        try:
            term = self.ontology[hpo_id]
            return [str(desc.id) for desc in term.subclasses() if str(desc.id) != hpo_id]
        except KeyError:
            return []
