import hail as hl
from typing import Callable, Dict, List, Optional

from hail.methods.misc import require_col_key_str, require_row_key_variant

from .concordance import (
    ConcordanceTable,
    global_concordance,
    sample_concordance,
    variant_concordance
)
from .utils import full_outer_join, JoinType


class JointCallSet:
    """A Python object representing the join between two MatrixTables.

    The JointCallSet supports two operations:

    - Concordance
    - Merge (Not Implemented Yet)

    Examples:

        # Necessary imports

        >>> from bge_toolkit.qc import ALL_AGG, ALL_GROUP, Statistic

        >>> exo = 'gs://my-bucket/exome_joint_called.mt'
        >>> imp = 'gs://my-bucket/imputed_with_glimpse.mt'
        >>> joint_callset = JointCallSet(exo, imp)

        # Define a custom binning function

        >>> def mac_bin(s: hl.StructExpression):
        ...     mac_bin = (hl.case()
        ...       .when(s.AC <= 1, '1')
        ...       .when(s.AC <= 5, '2-5')
        ...       .when(s.AC <= 10, '6-10')
        ...       .default('10+'))
        ...     return mac_bin

        # Import a binning function

        >>> from bge_toolkit.common import dp_bin

        # Specify aggregation variables

        >>> exo_row_group_by_var = {'exo_mac': mac_bin}
        >>> exo_entry_agg_fields = {'DP': dp_bin}
        >>> exo_col_group_by_var = {'pop': lambda s: s.pop}

        # Compute global concordance and make a plot

        >>> global_conc = joint_callset.global_concordance(exo_row_group_by_var=exo_row_group_by_var,
        ...                                                exo_col_group_by_var=exo_col_group_by_var,
        ...                                                exo_entry_agg_fields=exo_entry_agg_fields)
        >>> global_view = global_conc.get_view(group_names=['exo_mac', 'pop'], agg_names=['DP'], ordering={'exo_mac': ['1', '2-5', '6-10', '10+']})
        >>> global_view.plot(statistic=Statistic.NONREF_CONCORDANCE)
        >>> global_df = global_view.result()

        # Compute variant concordance and export results for a specific variant

        >>> variant_conc = joint_callset.variant_concordance(exo_col_group_by_var=exo_col_group_by_var,
        ...                                                  exo_entry_agg_fields=exo_entry_agg_fields)
        >>> variant_view = variant_conc.get_view(key=dict(locus=hl.Locus('chr1', 3452452, reference_genome='GRCh38'), alleles=['A', 'G']),
        ...                                      group_names=[ALL_GROUP, 'pop'],
        ...                                      agg_names=[ALL_AGG, 'DP'])
        >>> variant_view.export('gs://my-bucket/results.tsv')

        # Compute sample concordance and get the results for all samples

        >>> sample_conc = joint_callset.sample_concordance()
        >>> sample_view = sample_conc.get_view()
        >>> sample_view.result()

        # When done, close to unpersist any cached data

        >>> joint_callset.close()

    Args:
        exo (hl.MatrixTable): A MatrixTable that was generated from whole-exome sequencing.
        imp (hl.MatrixTable): A MatrixTable that was generated from imputation.

    """
    def __init__(self, exo: hl.MatrixTable, imp: hl.MatrixTable):
        require_col_key_str(exo, 'exo')
        require_col_key_str(imp, 'imp')

        require_row_key_variant(exo, 'exo')
        require_row_key_variant(imp, 'imp')

        self.exo = exo
        self.imp = imp

        self._sample_overlaps: Optional[hl.Table] = None
        self._variant_overlaps: Optional[hl.Table] = None
        self._conc_tables: List[ConcordanceTable] = []

    def close(self):
        """Cleanup any cached datasets before exiting."""
        if self._sample_overlaps is not None:
            self._sample_overlaps.unpersist()
        if self._variant_overlaps is not None:
            self._variant_overlaps.unpersist()
        for conc_t in self._conc_tables:
            conc_t.close()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __enter__(self):
        return self

    @property
    def n_sample_overlaps(self) -> int:
        """Get the number of overlapping samples between the exome and imputation datasets.

        Returns:
            int: The number of overlaps.
        """
        sample_overlaps = self.sample_overlaps()
        n_overlaps = sample_overlaps.aggregate(hl.struct(n_overlaps=hl.agg.count_where(sample_overlaps.EXOME & sample_overlaps.IMPUTATION)))
        return n_overlaps.n_overlaps

    def sample_overlaps(self) -> hl.Table:
        """Get a HailTable with the union of all samples along with their membership in each dataset.

        Returns:
            hl.Table: A table with all samples with two columns: EXOME and IMPUTATION. Both columns are booleans which represent a sample's member in each dataset.
        """
        if self._sample_overlaps is None:
            exo_samples = self.exo.cols().select()
            exo_samples = exo_samples.annotate(EXOME=True)

            imp_samples = self.imp.cols().select()
            imp_samples = imp_samples.annotate(IMPUTATION=True)

            result = exo_samples.join(imp_samples, how='outer')
            result = result.annotate(EXOME=hl.coalesce(result.EXOME, False),
                                     IMPUTATION=hl.coalesce(result.IMPUTATION, False))

            result = result.persist()

            self._sample_overlaps = result

        return self._sample_overlaps

    @property
    def n_variant_overlaps(self) -> int:
        """Get the number of overlapping variants between the exome and imputation datasets.

        Returns:
            int: The number of overlaps.
        """
        variant_overlaps = self.variant_overlaps()
        return variant_overlaps.count()

    def variant_overlaps(self) -> hl.Table:
        """Get a HailTable with all overlapping variants.

        Returns:
            hl.Table: A table with all overlapping variants.
        """
        if self._variant_overlaps is None:
            exo_variants = self.exo.rows().select()
            imp_variants = self.imp.rows().select()
            result = exo_variants.join(imp_variants, how='inner')
            result = result.persist()
            self._variant_overlaps = result
        return self._variant_overlaps

    def full_outer_join(self, join_type: JoinType = JoinType.outer) -> hl.MatrixTable:
        """Return a Hail MatrixTable with a full outer join.

        Args:
            join_type (JoinType): How to join the two datasets. One of "inner" or "outer".

        Returns:
            hl.MatrixTable: Equivalent of running hl.experimental.full_outer_join_mt()
        """
        return full_outer_join(self.exo, self.imp, join_type)

    def global_concordance(self,
                           *,
                           exo_row_group_by_var: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                           imp_row_group_by_var: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                           exo_col_group_by_var: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                           imp_col_group_by_var: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                           exo_row_agg_fields: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                           imp_row_agg_fields: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                           exo_col_agg_fields: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                           imp_col_agg_fields: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                           exo_entry_agg_fields: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                           imp_entry_agg_fields: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                           join_type: JoinType = JoinType.inner) -> ConcordanceTable:
        """Compute the global concordance with various grouping and aggregation variables.

        Examples:

            # Import shortcuts for ungrouped groupings and aggregations

            >>> from bge_toolkit.qc.utils import ALL_AGG, ALL_GROUP

            # Define a custom binning function

            >>> def mac_bin(s: hl.StructExpression):
            ...     mac_bin = (hl.case()
            ...       .when(s.AC <= 1, '1')
            ...       .when(s.AC <= 5, '2-5')
            ...       .when(s.AC <= 10, '6-10')
            ...       .default('10+'))
            ...     return mac_bin

            # Or import a binning function

            >>> from bge_toolkit.common import dp_bin

            # Specify aggregation variables

            >>> exo_row_group_by_var = {'exo_mac': mac_bin}
            >>> exo_entry_agg_fields = {'DP': dp_bin}
            >>> exo_col_group_by_var = {'pop': lambda s: s.pop}

            # Compute global concordance and make a plot

            >>> global_conc = joint_callset.global_concordance(exo_row_group_by_var=exo_row_group_by_var,
            ...                                                exo_col_group_by_var=exo_col_group_by_var,
            ...                                                exo_entry_agg_fields=exo_entry_agg_fields)
            >>> global_view = global_conc.get_view(group_names=['exo_mac', 'pop'], agg_names=['DP', ALL_AGG], ordering={'exo_mac': ['1', '2-5', '6-10', '10+']})
            >>> global_view.plot()

        Notes:
            The total ungrouped concordance can be accessed with ALL_GROUP and ALL_AGG.

        Args:
            exo_row_group_by_var (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.left_row (exome) in a full outer join and convert it to a binned variable for grouping.
            imp_row_group_by_var (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.right_row (imputation) in a full outer join and convert it to a binned variable for grouping.
            exo_col_group_by_var (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.left_col (exome) in a full outer join and convert it to a binned variable for grouping.
            imp_col_group_by_var (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.right_col (imputation) in a full outer join and convert it to a binned variable for grouping.
            exo_row_agg_fields (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.left_row (exome) in a full outer join and convert it to a binned variable for aggregation.
            imp_row_agg_fields (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.right_row (imputation) in a full outer join and convert it to a binned variable for aggregation.
            exo_col_agg_fields (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.left_col (exome) in a full outer join and convert it to a binned variable for aggregation.
            imp_col_agg_fields (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.right_col (imputation) in a full outer join and convert it to a binned variable for aggregation.
            exo_entry_agg_fields (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.left_entry (exome) in a full outer join and convert it to a binned variable for aggregation.
            imp_entry_agg_fields (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.right_entry (imputation) in a full outer join and convert it to a binned variable for aggregation.
            join_type (JoinType): Specify which variants to calculate concordance from. Can be one of "inner" or "outer". "outer" is substantially more expensive.

        Returns:
            ConcordanceTable: A table with the global concordance results.
        """
        conc_t = global_concordance(self.exo,
                                    self.imp,
                                    exo_row_group_by_var=exo_row_group_by_var or {},
                                    imp_row_group_by_var=imp_row_group_by_var or {},
                                    exo_col_group_by_var=exo_col_group_by_var or {},
                                    imp_col_group_by_var=imp_col_group_by_var or {},
                                    exo_row_agg_fields=exo_row_agg_fields or {},
                                    imp_row_agg_fields=imp_row_agg_fields or {},
                                    exo_col_agg_fields=exo_col_agg_fields or {},
                                    imp_col_agg_fields=imp_col_agg_fields or {},
                                    exo_entry_agg_fields=exo_entry_agg_fields or {},
                                    imp_entry_agg_fields=imp_entry_agg_fields or {},
                                    join_type=join_type)
        self._conc_tables.append(conc_t)
        return conc_t

    def variant_concordance(self,
                            *,
                            exo_col_group_by_var: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                            imp_col_group_by_var: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                            exo_col_agg_fields: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                            imp_col_agg_fields: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                            exo_entry_agg_fields: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                            imp_entry_agg_fields: Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]] = None,
                            join_type: JoinType = JoinType.inner) -> ConcordanceTable:
        """Compute the variant concordance with various grouping and aggregation variables.

        Examples:

            # Import shortcuts for ungrouped groupings and aggregations

            >>> from bge_toolkit.qc.utils import ALL_AGG, ALL_GROUP

            # Define a custom binning function

            >>> def dp_bin(s: hl.StructExpression):
            ...     x = (hl.case()
            ...     .when(s.DP <= 10, 10)
            ...     .when(s.DP <= 20, 20)
            ...     .when(s.DP <= 30, 30)
            ...     .when(s.DP <= 40, 40)
            ...     .when(s.DP <= 50, 50)
            ...     .when(s.DP <= 60, 60)
            ...     .when(s.DP <= 70, 70)
            ...     .when(s.DP <= 80, 80)
            ...     .when(s.DP <= 90, 90)
            ...     .default(100))
            ...     return x

            # Specify aggregation variables

            >>> exo_entry_agg_fields = {'DP': dp_bin}
            >>> exo_col_group_by_var = {'pop': lambda s: s.pop}

            # Compute variant concordance and export results for a specific variant

            >>> variant_conc = joint_callset.variant_concordance(exo_col_group_by_var=exo_col_group_by_var,
            ...                                                  exo_entry_agg_fields=exo_entry_agg_fields)
            >>> variant_view = variant_conc.get_view(key=dict(locus=hl.Locus('chr1', 3452452, reference_genome='GRCh38'), alleles=['A', 'G']),
            ...                                      group_names=[ALL_GROUP, 'pop'],
            ...                                      agg_names=[ALL_AGG, 'DP'])
            >>> variant_view.export('gs://my-bucket/results.tsv')

        Notes:
            The total ungrouped concordance can be accessed with ALL_GROUP and ALL_AGG.

        Args:
            exo_col_group_by_var (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.left_col (exome) in a full outer join and convert it to a binned variable for grouping.
            imp_col_group_by_var (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.right_col (imputation) in a full outer join and convert it to a binned variable for grouping.
            exo_col_agg_fields (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.left_col (exome) in a full outer join and convert it to a binned variable for aggregation.
            imp_col_agg_fields (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.right_col (imputation) in a full outer join and convert it to a binned variable for aggregation.
            exo_entry_agg_fields (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.left_entry (exome) in a full outer join and convert it to a binned variable for aggregation.
            imp_entry_agg_fields (Optional[Dict[str, Callable[[hl.StructExpression], hl.Expression]]]): An optional mapping from a grouping variable name to a
             function that specifies how to take a struct representing mt.right_entry (imputation) in a full outer join and convert it to a binned variable for aggregation.
            join_type (JoinType): Specify which variants to calculate concordance from. Can be one of "inner" or "outer".

        Returns:
            ConcordanceTable: A table with the concordance results for all variants.
        """
        conc_t = variant_concordance(self.exo,
                                     self.imp,
                                     exo_col_group_by_var=exo_col_group_by_var or {},
                                     imp_col_group_by_var=imp_col_group_by_var or {},
                                     exo_col_agg_fields=exo_col_agg_fields or {},
                                     imp_col_agg_fields=imp_col_agg_fields or {},
                                     exo_entry_agg_fields=exo_entry_agg_fields or {},
                                     imp_entry_agg_fields=imp_entry_agg_fields or {},
                                     join_type=join_type)
        self._conc_tables.append(conc_t)
        return conc_t

    def sample_concordance(self, *, join_type: JoinType = JoinType.inner) -> ConcordanceTable:
        """Compute the global sample concordance.

        Examples:

            # Import shortcuts for ungrouped groupings and aggregations

            >>> from bge_toolkit.qc.utils import ALL_AGG, ALL_GROUP


            # Compute sample concordance and make a plot for a specific sample

            >>> sample_conc = joint_callset.sample_concordance()
            >>> sample_view = sample_conc.get_view()
            >>> sample_view.result()

        Args:
            join_type (JoinType): Specify which variants to calculate concordance from. Can be one of "inner" or "outer".

        Returns:
            ConcordanceTable: A table with the concordance results for all samples.
        """
        conc_t = sample_concordance(self.exo, self.imp, join_type=join_type)
        self._conc_tables.append(conc_t)
        return conc_t
