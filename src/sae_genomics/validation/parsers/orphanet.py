"""Orphanet rare disease database parser."""

from pathlib import Path
from typing import Iterator, Optional
import xml.etree.ElementTree as ET

from .base import DatabaseParser, GeneDiseaseAssociation


class OrphanetParser(DatabaseParser):
    """Parser for Orphanet rare disease gene associations.

    Orphanet provides expert-curated information on rare diseases
    and their genetic associations.

    File format: XML with hierarchical structure containing:
        - Disorder information (OrphaCode, name)
        - Gene associations (gene symbol, association type, status)
    """

    def __init__(self, data_dir: Path):
        """Initialize Orphanet parser.

        Args:
            data_dir: Directory containing Orphanet XML files
        """
        super().__init__(data_dir)

        # Find XML file (filename may vary by release date)
        xml_files = list(self.data_dir.glob('en_product6*.xml'))

        if not xml_files:
            raise FileNotFoundError(
                f"Orphanet gene-disease file not found in {self.data_dir}. "
                f"Expected: en_product6*.xml"
            )

        self.xml_file = xml_files[0]

    def parse_gene_disease_associations(self) -> Iterator[GeneDiseaseAssociation]:
        """Parse Orphanet gene-rare disease associations.

        Yields:
            GeneDiseaseAssociation objects
        """
        # Parse XML
        tree = ET.parse(self.xml_file)
        root = tree.getroot()

        # Iterate through disorders
        for disorder in root.findall('.//Disorder'):
            # Extract disorder information
            orpha_code_elem = disorder.find('.//OrphaCode')
            if orpha_code_elem is None:
                continue

            orpha_code = orpha_code_elem.text
            disease_id = f"ORPHA:{orpha_code}"

            # Get disorder name
            name_elem = disorder.find('.//Name')
            disease_name = name_elem.text if name_elem is not None else orpha_code

            # Get gene associations
            gene_list = disorder.find('.//DisorderGeneAssociationList')
            if gene_list is None:
                continue

            for gene_assoc in gene_list.findall('.//DisorderGeneAssociation'):
                # Extract gene information
                gene = gene_assoc.find('.//Gene')
                if gene is None:
                    continue

                gene_symbol_elem = gene.find('.//Symbol')
                if gene_symbol_elem is None:
                    continue

                gene_symbol = gene_symbol_elem.text

                # Get association type (e.g., "Disease-causing germline mutation(s)")
                assoc_type_elem = gene_assoc.find('.//DisorderGeneAssociationType/Name')
                assoc_type = assoc_type_elem.text if assoc_type_elem is not None else None

                # Get association status (e.g., "Assessed", "Not yet assessed")
                assoc_status_elem = gene_assoc.find('.//DisorderGeneAssociationStatus/Name')
                assoc_status = assoc_status_elem.text if assoc_status_elem is not None else None

                # Only include validated associations
                if assoc_status and 'not' in assoc_status.lower():
                    continue

                # Map association type to evidence level
                evidence_level = self._map_association_type(assoc_type)

                yield GeneDiseaseAssociation(
                    gene_id=gene_symbol,  # Will be mapped to Ensembl later
                    gene_symbol=gene_symbol,
                    disease_id=disease_id,
                    disease_name=disease_name,
                    source='orphanet',
                    evidence_level=evidence_level,
                    metadata={
                        'association_type': assoc_type,
                        'association_status': assoc_status,
                    }
                )

    def _map_association_type(self, assoc_type: Optional[str]) -> str:
        """Map Orphanet association type to standardized evidence level.

        Args:
            assoc_type: Orphanet association type

        Returns:
            Standardized evidence level
        """
        if not assoc_type:
            return 'unknown'

        assoc_lower = assoc_type.lower()

        if 'disease-causing' in assoc_lower or 'pathogenic' in assoc_lower:
            return 'definitive'
        elif 'major susceptibility' in assoc_lower:
            return 'strong'
        elif 'modifying' in assoc_lower or 'biomarker' in assoc_lower:
            return 'moderate'
        elif 'candidate' in assoc_lower:
            return 'limited'
        else:
            return 'curated'

    def get_disorder_by_orpha_code(self, orpha_code: str) -> Optional[dict]:
        """Get disorder information by OrphaCode.

        Args:
            orpha_code: Orphanet disorder ID (e.g., "ORPHA:558" or "558")

        Returns:
            Dict with disorder information or None if not found
        """
        # Remove ORPHA: prefix if present
        if orpha_code.startswith('ORPHA:'):
            orpha_code = orpha_code[6:]

        tree = ET.parse(self.xml_file)
        root = tree.getroot()

        for disorder in root.findall('.//Disorder'):
            code_elem = disorder.find('.//OrphaCode')
            if code_elem is not None and code_elem.text == orpha_code:
                name_elem = disorder.find('.//Name')

                # Get associated genes
                genes = []
                gene_list = disorder.find('.//DisorderGeneAssociationList')
                if gene_list is not None:
                    for gene_assoc in gene_list.findall('.//DisorderGeneAssociation'):
                        gene = gene_assoc.find('.//Gene')
                        if gene is not None:
                            symbol_elem = gene.find('.//Symbol')
                            if symbol_elem is not None:
                                genes.append(symbol_elem.text)

                return {
                    'orpha_code': orpha_code,
                    'name': name_elem.text if name_elem is not None else None,
                    'genes': genes,
                }

        return None

    def get_genes_for_disorder(self, orpha_code: str) -> list:
        """Get all genes associated with a disorder.

        Args:
            orpha_code: Orphanet disorder ID

        Returns:
            List of gene symbols
        """
        disorder = self.get_disorder_by_orpha_code(orpha_code)
        return disorder['genes'] if disorder else []
