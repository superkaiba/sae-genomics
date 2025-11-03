"""ClinVar variant-disease database parser."""

from pathlib import Path
from typing import Iterator, Optional, Set, Dict, List
import gzip
import re

from .base import DatabaseParser, GeneDiseaseAssociation


class ClinVarVariant:
    """Represents a ClinVar variant with clinical significance."""

    def __init__(
        self,
        variant_id: str,
        gene_symbol: str,
        disease_name: str,
        clinical_significance: str,
        review_status: str,
        phenotype_ids: List[str] = None,
    ):
        """Initialize ClinVar variant.

        Args:
            variant_id: ClinVar variant ID (e.g., "RCV000000012")
            gene_symbol: Gene symbol
            disease_name: Disease/phenotype name
            clinical_significance: Clinical significance (e.g., "Pathogenic")
            review_status: Review status (e.g., "reviewed by expert panel")
            phenotype_ids: List of phenotype IDs (MedGen, OMIM, etc.)
        """
        self.variant_id = variant_id
        self.gene_symbol = gene_symbol
        self.disease_name = disease_name
        self.clinical_significance = clinical_significance
        self.review_status = review_status
        self.phenotype_ids = phenotype_ids or []

    def is_pathogenic(self) -> bool:
        """Check if variant is pathogenic or likely pathogenic."""
        sig_lower = self.clinical_significance.lower()
        return 'pathogenic' in sig_lower and 'benign' not in sig_lower


class ClinVarParser(DatabaseParser):
    """Parser for ClinVar variant-disease associations.

    ClinVar aggregates information about genomic variation and its
    relationship to human health.

    File format: VCF (Variant Call Format) with INFO field containing:
        - GENEINFO: Gene symbol and ID
        - CLNSIG: Clinical significance
        - CLNDN: Disease name
        - CLNDISDB: Disease database IDs
        - CLNREVSTAT: Review status
    """

    def __init__(
        self,
        data_dir: Path,
        pathogenic_only: bool = True,
        reviewed_only: bool = False,
    ):
        """Initialize ClinVar parser.

        Args:
            data_dir: Directory containing ClinVar VCF files
            pathogenic_only: Only include pathogenic/likely pathogenic variants
            reviewed_only: Only include variants with expert review
        """
        super().__init__(data_dir)
        self.pathogenic_only = pathogenic_only
        self.reviewed_only = reviewed_only

        # Find VCF file (may be gzipped)
        vcf_files = (
            list(self.data_dir.glob('clinvar*.vcf.gz')) +
            list(self.data_dir.glob('clinvar*.vcf'))
        )

        if not vcf_files:
            raise FileNotFoundError(
                f"ClinVar VCF file not found in {self.data_dir}. "
                f"Expected: clinvar*.vcf or clinvar*.vcf.gz"
            )

        self.vcf_file = vcf_files[0]
        self.is_gzipped = self.vcf_file.suffix == '.gz'

    def parse_variants(self) -> Iterator[ClinVarVariant]:
        """Parse ClinVar variants from VCF.

        Yields:
            ClinVarVariant objects
        """
        open_func = gzip.open if self.is_gzipped else open
        mode = 'rt' if self.is_gzipped else 'r'

        with open_func(self.vcf_file, mode) as f:
            for line in f:
                # Skip header lines
                if line.startswith('#'):
                    continue

                # Parse VCF line
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue

                info_field = parts[7]

                # Parse INFO field
                info = self._parse_info_field(info_field)

                # Extract gene information
                gene_info = info.get('GENEINFO', '')
                if not gene_info:
                    continue

                # GENEINFO format: "GENE1:ID1|GENE2:ID2"
                genes = [g.split(':')[0] for g in gene_info.split('|') if ':' in g]
                if not genes:
                    continue

                # Extract clinical significance
                clnsig = info.get('CLNSIG', '')
                if not clnsig:
                    continue

                # Filter by pathogenicity
                if self.pathogenic_only:
                    if not self._is_pathogenic(clnsig):
                        continue

                # Extract disease information
                disease_name = info.get('CLNDN', '')
                if not disease_name or disease_name == '.':
                    continue

                # Clean disease name (may have multiple separated by |)
                diseases = [d.strip() for d in disease_name.split('|') if d.strip()]

                # Extract phenotype IDs
                disease_db = info.get('CLNDISDB', '')
                phenotype_ids = self._parse_phenotype_ids(disease_db)

                # Extract review status
                review_status = info.get('CLNREVSTAT', '')

                # Filter by review status
                if self.reviewed_only:
                    if not self._has_expert_review(review_status):
                        continue

                # Extract variant ID
                variant_id = parts[2] if parts[2] != '.' else f"{parts[0]}:{parts[1]}"

                # Yield variant for each gene-disease pair
                for gene in genes:
                    for disease in diseases:
                        yield ClinVarVariant(
                            variant_id=variant_id,
                            gene_symbol=gene,
                            disease_name=disease,
                            clinical_significance=clnsig,
                            review_status=review_status,
                            phenotype_ids=phenotype_ids,
                        )

    def parse_gene_disease_associations(self) -> Iterator[GeneDiseaseAssociation]:
        """Parse ClinVar gene-disease associations.

        Aggregates variants by gene-disease pairs and counts pathogenic variants.

        Yields:
            GeneDiseaseAssociation objects
        """
        # Aggregate variants by gene-disease
        gene_disease_variants: Dict[tuple, List[ClinVarVariant]] = {}

        for variant in self.parse_variants():
            key = (variant.gene_symbol, variant.disease_name)
            if key not in gene_disease_variants:
                gene_disease_variants[key] = []
            gene_disease_variants[key].append(variant)

        # Create associations
        for (gene_symbol, disease_name), variants in gene_disease_variants.items():
            # Count pathogenic variants
            n_pathogenic = sum(1 for v in variants if v.is_pathogenic())

            # Get unique phenotype IDs
            phenotype_ids = set()
            for v in variants:
                phenotype_ids.update(v.phenotype_ids)

            # Determine evidence level based on review status
            review_statuses = [v.review_status for v in variants]
            evidence_level = self._aggregate_evidence_level(review_statuses)

            # Use first disease ID or disease name
            disease_id = list(phenotype_ids)[0] if phenotype_ids else disease_name

            yield GeneDiseaseAssociation(
                gene_id=gene_symbol,  # Will be mapped to Ensembl later
                gene_symbol=gene_symbol,
                disease_id=disease_id,
                disease_name=disease_name,
                source='clinvar',
                evidence_level=evidence_level,
                score=float(n_pathogenic),  # Number of pathogenic variants
                metadata={
                    'n_variants': len(variants),
                    'n_pathogenic': n_pathogenic,
                    'phenotype_ids': list(phenotype_ids),
                }
            )

    def _parse_info_field(self, info_field: str) -> Dict[str, str]:
        """Parse VCF INFO field into dictionary.

        Args:
            info_field: INFO field string

        Returns:
            Dict mapping INFO keys to values
        """
        info = {}
        for item in info_field.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                info[key] = value
            else:
                info[item] = True
        return info

    def _parse_phenotype_ids(self, disease_db: str) -> List[str]:
        """Parse phenotype IDs from CLNDISDB field.

        Args:
            disease_db: CLNDISDB field value (e.g., "MedGen:C0123456,OMIM:123456")

        Returns:
            List of phenotype IDs
        """
        if not disease_db or disease_db == '.':
            return []

        ids = []
        for entry in disease_db.split('|'):
            for db_entry in entry.split(','):
                if ':' in db_entry:
                    ids.append(db_entry.strip())
        return ids

    def _is_pathogenic(self, clnsig: str) -> bool:
        """Check if clinical significance indicates pathogenicity.

        Args:
            clnsig: Clinical significance string

        Returns:
            True if pathogenic or likely pathogenic
        """
        sig_lower = clnsig.lower()
        return 'pathogenic' in sig_lower and 'benign' not in sig_lower

    def _has_expert_review(self, review_status: str) -> bool:
        """Check if variant has expert review.

        Args:
            review_status: Review status string

        Returns:
            True if reviewed by expert panel
        """
        if not review_status:
            return False
        status_lower = review_status.lower()
        return 'expert' in status_lower or 'reviewed' in status_lower

    def _aggregate_evidence_level(self, review_statuses: List[str]) -> str:
        """Aggregate evidence level from multiple review statuses.

        Args:
            review_statuses: List of review status strings

        Returns:
            Aggregated evidence level
        """
        # Check for expert review
        for status in review_statuses:
            if self._has_expert_review(status):
                return 'definitive'

        # Check for multiple submitters
        for status in review_statuses:
            if 'multiple' in status.lower():
                return 'strong'

        # Default to curated
        return 'curated'

    def get_genes_with_pathogenic_variants(self, min_variants: int = 1) -> Set[str]:
        """Get genes with pathogenic variants.

        Args:
            min_variants: Minimum number of pathogenic variants

        Returns:
            Set of gene symbols
        """
        gene_counts = {}

        for variant in self.parse_variants():
            if variant.is_pathogenic():
                gene_counts[variant.gene_symbol] = gene_counts.get(variant.gene_symbol, 0) + 1

        return {gene for gene, count in gene_counts.items() if count >= min_variants}
