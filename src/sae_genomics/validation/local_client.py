"""Local validation client for querying downloaded databases."""

from pathlib import Path
from typing import Dict, List, Set, Optional
from collections import defaultdict

from rich.console import Console
from rich.progress import track

from .parsers import (
    GenCCParser,
    HPOParser,
    GWASCatalogParser,
    STRINGParser,
    OrphanetParser,
    ClinVarParser,
    MonarchParser,
    OpenTargetsParser,
)
from .enrichment import EnrichmentAnalyzer, EnrichmentResult


console = Console()


class LocalValidationClient:
    """Client for querying local validation databases.

    This class provides a unified interface for:
    1. Loading and parsing multiple databases
    2. Running enrichment analysis
    3. Querying gene-disease/phenotype associations
    """

    def __init__(
        self,
        databases_dir: Path,
        load_on_init: bool = True,
        verbose: bool = True,
    ):
        """Initialize local validation client.

        Args:
            databases_dir: Directory containing downloaded databases
            load_on_init: Whether to load databases on initialization
            verbose: Whether to print loading progress
        """
        self.databases_dir = Path(databases_dir)
        self.verbose = verbose

        # Database parsers
        self.parsers: Dict = {}

        # Parsed data structures
        self.data: Dict = {
            'gencc': {'terms_to_genes': {}, 'term_names': {}},
            'hpo_phenotypes': {'terms_to_genes': {}, 'term_names': {}},
            'hpo_diseases': {'terms_to_genes': {}, 'term_names': {}},
            'gwas': {'terms_to_genes': {}, 'term_names': {}},
            'orphanet': {'terms_to_genes': {}, 'term_names': {}},
            'clinvar': {'terms_to_genes': {}, 'term_names': {}},
            'monarch': {'terms_to_genes': {}, 'term_names': {}},
            'open_targets': {'terms_to_genes': {}, 'term_names': {}},
            # NOTE: all_genes uses gene_symbol (most portable across databases)
            # For Ensembl ID-based enrichment, gene symbols must be mapped to Ensembl IDs
            # TODO: Add symbol-to-Ensembl mapping infrastructure for proper ID-based enrichment
            'all_genes': set(),
        }

        # STRING is special (protein interactions, not gene-disease)
        self.string_parser: Optional[STRINGParser] = None

        if load_on_init:
            self.load_databases()

    def _print(self, message: str, style: str = ""):
        """Print message if verbose mode is on."""
        if self.verbose:
            if style:
                console.print(f"[{style}]{message}[/{style}]")
            else:
                console.print(message)

    def load_databases(self):
        """Load and parse all available databases."""
        self._print("\n[bold]Loading validation databases...[/bold]\n")

        # Load GenCC
        gencc_dir = self.databases_dir / 'gencc'
        if gencc_dir.exists():
            self._load_gencc(gencc_dir)

        # Load HPO
        hpo_dir = self.databases_dir / 'hpo'
        if hpo_dir.exists():
            self._load_hpo(hpo_dir)

        # Load GWAS Catalog
        gwas_dir = self.databases_dir / 'gwas_catalog'
        if gwas_dir.exists():
            self._load_gwas(gwas_dir)

        # Load Orphanet
        orphanet_dir = self.databases_dir / 'orphanet'
        if orphanet_dir.exists():
            self._load_orphanet(orphanet_dir)

        # Load STRING (protein interactions)
        string_dir = self.databases_dir / 'string'
        if string_dir.exists():
            self._load_string(string_dir)

        # Load ClinVar
        clinvar_dir = self.databases_dir / 'clinvar'
        if clinvar_dir.exists():
            self._load_clinvar(clinvar_dir)

        # Load Monarch
        monarch_dir = self.databases_dir / 'monarch'
        if monarch_dir.exists():
            self._load_monarch(monarch_dir)

        # Load Open Targets
        open_targets_dir = self.databases_dir / 'open_targets'
        if open_targets_dir.exists():
            self._load_open_targets(open_targets_dir)

        # Summary
        n_genes = len(self.data['all_genes'])
        self._print(f"\n[green]✓ Loaded {n_genes} unique genes across all databases[/green]\n")

    def _load_gencc(self, gencc_dir: Path):
        """Load GenCC database."""
        self._print("[blue]Loading GenCC...[/blue]")
        parser = GenCCParser(gencc_dir)
        self.parsers['gencc'] = parser

        terms = defaultdict(set)
        names = {}

        for assoc in track(
            list(parser.parse_gene_disease_associations()),
            description="Parsing GenCC",
            disable=not self.verbose,
        ):
            terms[assoc.disease_id].add(assoc.gene_symbol)
            names[assoc.disease_id] = assoc.disease_name
            self.data['all_genes'].add(assoc.gene_symbol)

        self.data['gencc']['terms_to_genes'] = dict(terms)
        self.data['gencc']['term_names'] = names

        self._print(f"[green]✓ Loaded {len(terms)} diseases from GenCC[/green]")

    def _load_hpo(self, hpo_dir: Path):
        """Load HPO database."""
        self._print("[blue]Loading HPO...[/blue]")
        parser = HPOParser(hpo_dir)
        self.parsers['hpo'] = parser

        # Phenotype associations
        pheno_terms = defaultdict(set)
        pheno_names = {}

        for assoc in track(
            list(parser.parse_gene_phenotype_associations()),
            description="Parsing HPO phenotypes",
            disable=not self.verbose,
        ):
            pheno_terms[assoc.phenotype_id].add(assoc.gene_symbol)
            pheno_names[assoc.phenotype_id] = assoc.phenotype_name
            self.data['all_genes'].add(assoc.gene_symbol)

        self.data['hpo_phenotypes']['terms_to_genes'] = dict(pheno_terms)
        self.data['hpo_phenotypes']['term_names'] = pheno_names

        # Disease associations
        disease_terms = defaultdict(set)
        disease_names = {}

        for assoc in track(
            list(parser.parse_gene_disease_associations()),
            description="Parsing HPO diseases",
            disable=not self.verbose,
        ):
            disease_terms[assoc.disease_id].add(assoc.gene_symbol)
            disease_names[assoc.disease_id] = assoc.disease_name
            self.data['all_genes'].add(assoc.gene_symbol)

        self.data['hpo_diseases']['terms_to_genes'] = dict(disease_terms)
        self.data['hpo_diseases']['term_names'] = disease_names

        self._print(
            f"[green]✓ Loaded {len(pheno_terms)} phenotypes and "
            f"{len(disease_terms)} diseases from HPO[/green]"
        )

    def _load_gwas(self, gwas_dir: Path):
        """Load GWAS Catalog."""
        self._print("[blue]Loading GWAS Catalog...[/blue]")
        parser = GWASCatalogParser(gwas_dir, p_value_threshold=5e-8)
        self.parsers['gwas'] = parser

        terms = defaultdict(set)
        names = {}

        for assoc in track(
            list(parser.parse_gene_disease_associations()),
            description="Parsing GWAS",
            disable=not self.verbose,
        ):
            terms[assoc.disease_id].add(assoc.gene_symbol)
            names[assoc.disease_id] = assoc.disease_name
            self.data['all_genes'].add(assoc.gene_symbol)

        self.data['gwas']['terms_to_genes'] = dict(terms)
        self.data['gwas']['term_names'] = names

        self._print(f"[green]✓ Loaded {len(terms)} traits from GWAS Catalog[/green]")


    def _load_orphanet(self, orphanet_dir: Path):
        """Load Orphanet database."""
        self._print("[blue]Loading Orphanet...[/blue]")
        parser = OrphanetParser(orphanet_dir)
        self.parsers['orphanet'] = parser

        terms = defaultdict(set)
        names = {}

        for assoc in track(
            list(parser.parse_gene_disease_associations()),
            description="Parsing Orphanet",
            disable=not self.verbose,
        ):
            terms[assoc.disease_id].add(assoc.gene_symbol)
            names[assoc.disease_id] = assoc.disease_name
            self.data['all_genes'].add(assoc.gene_symbol)

        self.data['orphanet']['terms_to_genes'] = dict(terms)
        self.data['orphanet']['term_names'] = names

        self._print(f"[green]✓ Loaded {len(terms)} rare diseases from Orphanet[/green]")


    def _load_string(self, string_dir: Path):
        """Load STRING database."""
        self._print("[blue]Loading STRING...[/blue]")
        self.string_parser = STRINGParser(string_dir, score_threshold=400)
        self._print("[green]✓ STRING database ready[/green]")


    def _load_clinvar(self, clinvar_dir: Path):
        """Load ClinVar database."""
        self._print("[blue]Loading ClinVar...[/blue]")
        parser = ClinVarParser(clinvar_dir, pathogenic_only=True)
        self.parsers['clinvar'] = parser

        terms = defaultdict(set)
        names = {}

        # ClinVar is large - limit parsing to avoid memory issues
        count = 0
        max_variants = 100000  # Limit to first 100k variants

        for assoc in parser.parse_gene_disease_associations():
            if count >= max_variants:
                break

            terms[assoc.disease_id].add(assoc.gene_symbol)
            names[assoc.disease_id] = assoc.disease_name
            self.data['all_genes'].add(assoc.gene_symbol)
            count += 1

            # Progress update every 10k
            if count % 10000 == 0 and self.verbose:
                console.print(f"  Parsed {count:,} variants...", end='\r')

        self.data['clinvar']['terms_to_genes'] = dict(terms)
        self.data['clinvar']['term_names'] = names

        self._print(f"[green]✓ Loaded {len(terms)} conditions from ClinVar ({count:,} variants)[/green]")


    def _load_monarch(self, monarch_dir: Path):
        """Load Monarch database."""
        self._print("[blue]Loading Monarch...[/blue]")
        parser = MonarchParser(monarch_dir, format='duckdb')
        self.parsers['monarch'] = parser

        # Monarch is large (1M+ nodes, 14M+ edges)
        # We'll load it on-demand during validation using build_gene_disease_index()
        # This method is called during validate_gene_set to query diseases for specific genes

        self._print("[green]✓ Monarch database ready (on-demand queries)[/green]")


    def _load_open_targets(self, open_targets_dir: Path):
        """Load Open Targets database."""
        self._print("[blue]Loading Open Targets...[/blue]")
        parser = OpenTargetsParser(open_targets_dir, score_threshold=0.2)
        self.parsers['open_targets'] = parser

        # Open Targets is large (39K diseases, 78K targets, millions of associations)
        # We'll load it on-demand during validation to avoid memory issues

        self._print("[green]✓ Open Targets database ready (on-demand queries)[/green]")


    def validate_gene_set(
        self,
        gene_set: Set[str],
        databases: Optional[List[str]] = None,
        significance_threshold: float = 0.05,
        min_overlap: int = 2,
    ) -> Dict[str, List[EnrichmentResult]]:
        """Run enrichment analysis for a gene set across databases.

        Args:
            gene_set: Set of gene symbols to validate
            databases: List of databases to query (default: all available)
            significance_threshold: FDR significance threshold
            min_overlap: Minimum overlap with term

        Returns:
            Dict mapping database names to lists of significant enrichment results
        """
        if databases is None:
            databases = ['gencc', 'hpo_phenotypes', 'hpo_diseases', 'gwas', 'orphanet', 'clinvar', 'monarch', 'open_targets']

        # Filter to available databases
        available_dbs = []
        for db in databases:
            if db in ['monarch', 'open_targets']:
                # Monarch and Open Targets are query-based, check if parser exists
                if db in self.parsers:
                    available_dbs.append(db)
            elif db in self.data and self.data[db]['terms_to_genes']:
                available_dbs.append(db)

        databases = available_dbs

        # Initialize enrichment analyzer
        analyzer = EnrichmentAnalyzer(
            background_genes=self.data['all_genes'],
            significance_threshold=significance_threshold,
            min_overlap=min_overlap,
        )

        # Run enrichment for each database
        results = {}

        for db_name in databases:
            # Handle Monarch specially (query-based)
            if db_name == 'monarch':
                if 'monarch' in self.parsers:
                    # Query diseases for all genes in set
                    associations = list(self.parsers['monarch'].query_genes_batch(gene_set))

                    # Build disease→genes index from associations
                    terms_to_genes = defaultdict(set)
                    term_names = {}

                    for assoc in associations:
                        terms_to_genes[assoc.disease_id].add(assoc.gene_symbol)
                        term_names[assoc.disease_id] = assoc.disease_name

                    if terms_to_genes:
                        enrichment_results = analyzer.test_multiple_terms(
                            query_genes=gene_set,
                            term_to_genes=dict(terms_to_genes),
                            term_to_name=term_names,
                            source='monarch',
                        )
                        significant = analyzer.filter_significant(enrichment_results)
                        results['monarch'] = significant
                continue

            # Handle Open Targets specially (query-based)
            if db_name == 'open_targets':
                if 'open_targets' in self.parsers:
                    # Query diseases for all genes in set
                    associations = list(self.parsers['open_targets'].query_genes_batch(gene_set))

                    # Build disease→genes index from associations
                    terms_to_genes = defaultdict(set)
                    term_names = {}

                    for assoc in associations:
                        terms_to_genes[assoc.disease_id].add(assoc.gene_symbol)
                        term_names[assoc.disease_id] = assoc.disease_name

                    if terms_to_genes:
                        enrichment_results = analyzer.test_multiple_terms(
                            query_genes=gene_set,
                            term_to_genes=dict(terms_to_genes),
                            term_to_name=term_names,
                            source='open_targets',
                        )
                        significant = analyzer.filter_significant(enrichment_results)
                        results['open_targets'] = significant
                continue

            db_data = self.data[db_name]

            if not db_data['terms_to_genes']:
                continue

            # Run enrichment
            enrichment_results = analyzer.test_multiple_terms(
                query_genes=gene_set,
                term_to_genes=db_data['terms_to_genes'],
                term_to_name=db_data['term_names'],
                source=db_name,
            )

            # Filter significant
            significant = analyzer.filter_significant(enrichment_results)
            results[db_name] = significant

        return results

    def get_protein_interactions(
        self,
        gene_set: Set[str],
        min_score: int = 400,
    ) -> Dict[str, Set[str]]:
        """Get protein-protein interaction network for a gene set.

        Args:
            gene_set: Set of gene symbols
            min_score: Minimum STRING confidence score

        Returns:
            Dict mapping each gene to its interactors within the gene set
        """
        if self.string_parser is None:
            return {}

        return self.string_parser.build_network(gene_set, min_score=min_score)

    def get_genes_for_disease(self, disease_name: str, database: str = 'gencc') -> Set[str]:
        """Get genes associated with a disease.

        Args:
            disease_name: Disease name or ID
            database: Database to query

        Returns:
            Set of gene symbols
        """
        if database not in self.data:
            return set()

        db_data = self.data[database]
        genes = set()

        # Search by ID
        if disease_name in db_data['terms_to_genes']:
            genes.update(db_data['terms_to_genes'][disease_name])

        # Search by name (case-insensitive)
        disease_lower = disease_name.lower()
        for term_id, term_name in db_data['term_names'].items():
            if disease_lower in term_name.lower():
                genes.update(db_data['terms_to_genes'][term_id])

        return genes

    def get_database_stats(self) -> Dict:
        """Get statistics about loaded databases.

        Returns:
            Dict with database statistics
        """
        stats = {
            'total_genes': len(self.data['all_genes']),
            'databases': {},
        }

        for db_name, db_data in self.data.items():
            if db_name == 'all_genes':
                continue

            if 'terms_to_genes' in db_data:
                stats['databases'][db_name] = {
                    'n_terms': len(db_data['terms_to_genes']),
                    'n_genes': len({
                        gene
                        for genes in db_data['terms_to_genes'].values()
                        for gene in genes
                    }),
                }

        return stats
