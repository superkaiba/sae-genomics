"""GenCC database parser."""

from pathlib import Path
from typing import Iterator
import csv

from .base import DatabaseParser, GeneDiseaseAssociation


class GenCCParser(DatabaseParser):
    """Parser for GenCC gene-disease validity classifications.

    GenCC provides expert-curated gene-disease validity classifications
    from multiple sources (ClinGen, Genomics England, etc.).

    File format: CSV with columns including:
        - gene_curie: Gene ID (HGNC:12345)
        - gene_symbol: Gene symbol
        - disease_curie: Disease ID (MONDO:0012345, Orphanet:123, etc.)
        - disease_title: Disease name
        - classification_title: Validity classification
        - moi_title: Mode of inheritance
        - submitter_title: Submitting organization
    """

    def __init__(self, data_dir: Path):
        """Initialize GenCC parser.

        Args:
            data_dir: Directory containing gencc-submissions.csv
        """
        super().__init__(data_dir)
        self.csv_file = self.data_dir / 'gencc-submissions.csv'

        if not self.csv_file.exists():
            raise FileNotFoundError(f"GenCC file not found: {self.csv_file}")

    def parse_gene_disease_associations(self) -> Iterator[GeneDiseaseAssociation]:
        """Parse GenCC gene-disease associations.

        Yields:
            GeneDiseaseAssociation objects
        """
        with open(self.csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)

            for row in reader:
                # Extract gene information
                gene_curie = row.get('gene_curie', '')  # e.g., "HGNC:12345"
                gene_symbol = row.get('gene_symbol', '')

                # Skip if missing essential info
                if not gene_symbol:
                    continue

                # Convert HGNC ID to Ensembl (will be done by gene ID mapper later)
                # For now, use gene symbol as identifier
                gene_id = gene_symbol  # Placeholder, will be mapped later

                # Extract disease information
                disease_curie = row.get('disease_curie', '')  # e.g., "MONDO:0012345"
                disease_title = row.get('disease_title', '')

                if not disease_title:
                    continue

                # Extract classification (evidence level)
                classification = row.get('classification_title', '')
                # Maps to: Definitive, Strong, Moderate, Limited, Disputed, Refuted

                # Mode of inheritance
                moi = row.get('moi_title', '')

                # Submitting organization
                submitter = row.get('submitter_title', '')

                # PMIDs if available
                pmids = []
                if 'pmid' in row and row['pmid']:
                    # PMIDs might be pipe-separated
                    pmids = [p.strip() for p in row['pmid'].split('|') if p.strip()]

                yield GeneDiseaseAssociation(
                    gene_id=gene_id,
                    gene_symbol=gene_symbol,
                    disease_id=disease_curie if disease_curie else disease_title,
                    disease_name=disease_title,
                    source='gencc',
                    evidence_level=classification.lower() if classification else None,
                    pmids=pmids,
                    metadata={
                        'mode_of_inheritance': moi,
                        'submitter': submitter,
                        'gene_curie': gene_curie,
                    }
                )
