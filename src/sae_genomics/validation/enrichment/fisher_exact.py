"""Fisher's exact test enrichment analysis with FDR correction."""

from typing import List, Dict, Set, Optional
from dataclasses import dataclass
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


@dataclass
class EnrichmentResult:
    """Result of enrichment analysis for a single term."""

    term_id: str  # Disease/phenotype ID
    term_name: str  # Human-readable name
    p_value: float  # Raw p-value from Fisher's exact test
    p_value_adjusted: float  # FDR-corrected p-value
    odds_ratio: float  # Odds ratio
    genes_in_query: int  # Number of query genes
    genes_in_term: int  # Number of genes associated with this term
    genes_overlap: int  # Number of genes in both query and term
    genes_overlap_list: List[str]  # List of overlapping gene symbols
    background_size: int  # Total number of genes in background
    source: str  # Database source

    @property
    def is_significant(self) -> bool:
        """Check if result is significant at alpha=0.05."""
        return self.p_value_adjusted < 0.05

    def to_dict(self) -> Dict:
        """Convert to dictionary with JSON-serializable types."""
        return {
            'term_id': self.term_id,
            'term_name': self.term_name,
            'p_value': float(self.p_value),
            'p_value_adjusted': float(self.p_value_adjusted),
            'odds_ratio': float(self.odds_ratio),
            'genes_in_query': int(self.genes_in_query),
            'genes_in_term': int(self.genes_in_term),
            'genes_overlap': int(self.genes_overlap),
            'genes_overlap_list': self.genes_overlap_list,
            'background_size': int(self.background_size),
            'source': self.source,
            'significant': bool(self.is_significant),
        }


class EnrichmentAnalyzer:
    """Perform statistical enrichment analysis using Fisher's exact test.

    This class performs overrepresentation analysis to identify
    diseases/phenotypes/pathways that are significantly associated
    with a query gene set.
    """

    def __init__(
        self,
        background_genes: Set[str],
        significance_threshold: float = 0.05,
        min_overlap: int = 2,
    ):
        """Initialize enrichment analyzer.

        Args:
            background_genes: Set of all genes in the background (universe)
            significance_threshold: Alpha level for significance (default: 0.05)
            min_overlap: Minimum number of overlapping genes required (default: 2)
        """
        self.background_genes = background_genes
        self.background_size = len(background_genes)
        self.significance_threshold = significance_threshold
        self.min_overlap = min_overlap

    def test_single_term(
        self,
        query_genes: Set[str],
        term_genes: Set[str],
        term_id: str,
        term_name: str,
        source: str = "unknown",
    ) -> Optional[EnrichmentResult]:
        """Test enrichment for a single term using Fisher's exact test.

        Constructs 2x2 contingency table:
                        | In Term | Not in Term |
        In Query        |    a    |      b      |
        Not in Query    |    c    |      d      |

        Args:
            query_genes: Set of query gene symbols
            term_genes: Set of gene symbols associated with the term
            term_id: Term identifier
            term_name: Human-readable term name
            source: Database source

        Returns:
            EnrichmentResult or None if insufficient overlap
        """
        # Calculate overlap
        overlap_genes = query_genes & term_genes & self.background_genes
        overlap_count = len(overlap_genes)

        # Skip if insufficient overlap
        if overlap_count < self.min_overlap:
            return None

        # Construct contingency table
        a = overlap_count  # In query AND in term
        b = len(query_genes & self.background_genes) - a  # In query, NOT in term
        c = len(term_genes & self.background_genes) - a  # NOT in query, but in term
        d = self.background_size - a - b - c  # NOT in query, NOT in term

        # Ensure no negative values (can happen with inconsistent gene sets)
        if b < 0 or c < 0 or d < 0:
            return None

        # Fisher's exact test (one-sided: test for over-representation)
        try:
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
        except Exception:
            return None

        return EnrichmentResult(
            term_id=term_id,
            term_name=term_name,
            p_value=p_value,
            p_value_adjusted=p_value,  # Will be corrected later
            odds_ratio=odds_ratio if odds_ratio != float('inf') else 999.0,
            genes_in_query=len(query_genes & self.background_genes),
            genes_in_term=len(term_genes & self.background_genes),
            genes_overlap=overlap_count,
            genes_overlap_list=sorted(list(overlap_genes)),
            background_size=self.background_size,
            source=source,
        )

    def test_multiple_terms(
        self,
        query_genes: Set[str],
        term_to_genes: Dict[str, Set[str]],
        term_to_name: Dict[str, str],
        source: str = "unknown",
        correction_method: str = 'fdr_bh',
    ) -> List[EnrichmentResult]:
        """Test enrichment for multiple terms with multiple testing correction.

        Args:
            query_genes: Set of query gene symbols
            term_to_genes: Dict mapping term IDs to sets of gene symbols
            term_to_name: Dict mapping term IDs to human-readable names
            source: Database source
            correction_method: Method for multiple testing correction
                ('fdr_bh' for Benjamini-Hochberg FDR, 'bonferroni', etc.)

        Returns:
            List of EnrichmentResult objects, sorted by adjusted p-value
        """
        # Test each term
        results = []
        for term_id, term_genes in term_to_genes.items():
            term_name = term_to_name.get(term_id, term_id)

            result = self.test_single_term(
                query_genes=query_genes,
                term_genes=term_genes,
                term_id=term_id,
                term_name=term_name,
                source=source,
            )

            if result is not None:
                results.append(result)

        # Apply multiple testing correction
        if len(results) > 0:
            p_values = [r.p_value for r in results]

            try:
                reject, p_values_adj, _, _ = multipletests(
                    p_values,
                    alpha=self.significance_threshold,
                    method=correction_method
                )

                # Update adjusted p-values
                for i, result in enumerate(results):
                    result.p_value_adjusted = p_values_adj[i]

            except Exception as e:
                import sys
                print(f"Warning: Multiple testing correction failed: {e}", file=sys.stderr)
                # If correction fails, use raw p-values
                for result in results:
                    result.p_value_adjusted = result.p_value

        # Sort by adjusted p-value
        results.sort(key=lambda x: x.p_value_adjusted)

        return results

    def filter_significant(
        self,
        results: List[EnrichmentResult],
        threshold: Optional[float] = None,
    ) -> List[EnrichmentResult]:
        """Filter results to only significant terms.

        Args:
            results: List of EnrichmentResult objects
            threshold: Significance threshold (default: use instance threshold)

        Returns:
            Filtered list of significant results
        """
        threshold = threshold or self.significance_threshold
        return [r for r in results if r.p_value_adjusted < threshold]
