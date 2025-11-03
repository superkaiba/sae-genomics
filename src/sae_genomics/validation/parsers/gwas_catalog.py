"""GWAS Catalog database parser."""

from pathlib import Path
from typing import Iterator, List, Set
import csv

from .base import DatabaseParser, GeneDiseaseAssociation


class GWASCatalogParser(DatabaseParser):
    """Parser for GWAS Catalog associations.

    GWAS Catalog provides genome-wide association study results linking
    SNPs to traits/diseases with associated genes.

    File format: TSV with columns including:
        - MAPPED_GENE: Gene symbols (comma-separated)
        - DISEASE/TRAIT: Phenotype/disease name
        - MAPPED_TRAIT_URI: Ontology term (EFO)
        - P-VALUE: Association p-value
        - PUBMEDID: Publication ID
    """

    def __init__(self, data_dir: Path, p_value_threshold: float = 5e-8):
        """Initialize GWAS Catalog parser.

        Args:
            data_dir: Directory containing GWAS Catalog files
            p_value_threshold: Significance threshold (default: 5e-8, genome-wide)
        """
        super().__init__(data_dir)
        self.p_value_threshold = p_value_threshold

        # Try both possible filenames
        possible_files = [
            self.data_dir / 'gwas-catalog-associations_ontology-annotated.tsv',
            self.data_dir / 'gwas-catalog-associations.tsv',
        ]

        self.tsv_file = None
        for f in possible_files:
            if f.exists():
                self.tsv_file = f
                break

        if self.tsv_file is None:
            raise FileNotFoundError(
                f"GWAS Catalog file not found in {self.data_dir}. "
                f"Expected one of: {[f.name for f in possible_files]}"
            )

    def parse_gene_disease_associations(self) -> Iterator[GeneDiseaseAssociation]:
        """Parse GWAS Catalog gene-trait associations.

        Yields:
            GeneDiseaseAssociation objects for genome-wide significant associations
        """
        with open(self.tsv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')

            for row in reader:
                # Extract p-value
                p_value_str = row.get('P-VALUE', '')
                try:
                    p_value = float(p_value_str) if p_value_str else None
                except (ValueError, TypeError):
                    continue

                # Filter by significance threshold
                if p_value is None or p_value > self.p_value_threshold:
                    continue

                # Extract genes (can be comma-separated)
                genes_str = row.get('MAPPED_GENE', '')
                if not genes_str or genes_str in ['NR', 'intergenic']:
                    continue

                # Split multiple genes
                genes = [g.strip() for g in genes_str.split(',') if g.strip()]

                # Extract trait information
                trait_name = row.get('DISEASE/TRAIT', '')
                if not trait_name:
                    continue

                # Extract ontology term (EFO)
                trait_uri = row.get('MAPPED_TRAIT_URI', '')
                trait_id = trait_uri.split('/')[-1] if trait_uri else trait_name

                # Extract study information
                pubmed_id = row.get('PUBMEDID', '')
                pmids = [pubmed_id] if pubmed_id else []

                # Extract risk allele and SNP information
                snp_id = row.get('SNPS', '')
                risk_allele = row.get('STRONGEST SNP-RISK ALLELE', '')

                # OR/Beta coefficient
                or_value_str = row.get('OR or BETA', '')
                try:
                    or_value = float(or_value_str) if or_value_str else None
                except (ValueError, TypeError):
                    or_value = None

                # Create association for each gene
                for gene_symbol in genes:
                    yield GeneDiseaseAssociation(
                        gene_id=gene_symbol,  # Will be mapped to Ensembl later
                        gene_symbol=gene_symbol,
                        disease_id=trait_id,
                        disease_name=trait_name,
                        source='gwas_catalog',
                        evidence_level='genome-wide significant',
                        score=p_value,  # Use p-value as score
                        pmids=pmids,
                        metadata={
                            'p_value': p_value,
                            'or_beta': or_value,
                            'snp_id': snp_id,
                            'risk_allele': risk_allele,
                            'trait_uri': trait_uri,
                        }
                    )

    def get_genes_for_trait(self, trait_name: str) -> Set[str]:
        """Get all genes associated with a specific trait.

        Args:
            trait_name: Trait/disease name or ID

        Returns:
            Set of gene symbols
        """
        genes = set()
        for assoc in self.parse_gene_disease_associations():
            if (trait_name.lower() in assoc.disease_name.lower() or
                trait_name in assoc.disease_id):
                genes.add(assoc.gene_symbol)
        return genes

    def get_traits_for_gene(self, gene_symbol: str) -> Set[str]:
        """Get all traits associated with a specific gene.

        Args:
            gene_symbol: Gene symbol

        Returns:
            Set of trait names
        """
        traits = set()
        for assoc in self.parse_gene_disease_associations():
            if assoc.gene_symbol.upper() == gene_symbol.upper():
                traits.add(assoc.disease_name)
        return traits
