from enum import Enum
import os
from typing import Any, Callable, Dict, List, Optional, Tuple
from tempfile import NamedTemporaryFile

from plotnine import (ggplot, aes, guides, guide_legend, geom_point, geom_line, scale_size_continuous,
                      scale_color_manual, theme_minimal, labs, theme, facet_grid, element_text, coord_fixed,
                      scale_color_brewer, element_rect, scale_color_discrete)

import pandas as pd

import hail as hl
from hail import ArrayExpression
import  hailtop.fs as hfs

from .utils import ALL_AGG, ALL_GROUP, JoinType, full_outer_join


class GenotypeType(Enum):
    MISSING = 0
    NO_CALL = 1
    HOM_REF = 2
    HET = 3
    HOM_ALT = 4


class Statistic(Enum):
    """
    The statistic to use.

    Attributes:
        nonref_concordance_rate: Non-reference concordance rate.
        concordance_rate: Concordance rate.
        f1_score: F1 score.
    """
    NONREF_CONCORDANCE = 'nonref_concordance_rate'
    CONCORDANCE = 'concordance_rate'
    F1_SCORE = 'f1_score'


def get_cast_type(typ: hl.HailType):
    if typ in (hl.tint64, hl.tint32):
        return int
    if typ in (hl.tfloat64, hl.tfloat32):
        return float
    if typ == hl.tbool:
        return bool
    return str


class ConcordanceArray:
    def __init__(self, data: hl.ArrayExpression):
        self._data = data

    def get_count(self, exo_gt: GenotypeType, imp_gt: GenotypeType):
        return self._data[exo_gt.value][imp_gt.value]

    def concordance_rate(self) -> hl.Float64Expression:
        n_concordant = self.n_concordant()
        total_obs = hl.sum(hl.map(lambda x: hl.sum(x[1:]), self._data[2:]))
        return hl.or_missing(total_obs > 0, n_concordant / total_obs)

    def nonref_concordance(self) -> hl.Float64Expression:
        on_diag = (self.get_count(GenotypeType.HET, GenotypeType.HET) +
                   self.get_count(GenotypeType.HOM_ALT, GenotypeType.HOM_ALT))

        total_obs = hl.sum(hl.map(lambda x: hl.sum(x[1:]), self._data[3:]))

        return hl.or_missing(total_obs > 0, on_diag / total_obs)

    def true_pos(self) -> hl.Int32Expression:
        return self.get_count(GenotypeType.HET, GenotypeType.HET) + self.get_count(GenotypeType.HOM_ALT, GenotypeType.HOM_ALT)

    def true_neg(self) -> hl.Int32Expression:
        return self.get_count(GenotypeType.HOM_REF, GenotypeType.HOM_REF)

    def false_pos(self) -> hl.Int32Expression:
        return (self.get_count(GenotypeType.HOM_REF, GenotypeType.HET) +
                    self.get_count(GenotypeType.HOM_REF, GenotypeType.HOM_ALT) +
                    self.get_count(GenotypeType.HET, GenotypeType.HOM_ALT))

    def false_neg(self) -> hl.Int32Expression:
        return (self.get_count(GenotypeType.HET, GenotypeType.HOM_REF) +
                     self.get_count(GenotypeType.HOM_ALT, GenotypeType.HOM_REF) +
                     self.get_count(GenotypeType.HOM_ALT, GenotypeType.HET))

    def f1_score(self) -> hl.Float64Expression:
        true_pos = self.true_pos()
        false_pos = self.false_pos()
        false_neg = self.false_neg()

        prec = hl.or_missing(true_pos + false_pos > 0, true_pos / (true_pos + false_pos))
        recall = hl.or_missing(true_pos + false_neg > 0, true_pos / (true_pos + false_neg))

        return 2 * prec * recall / (prec + recall)

    def n_concordant(self) -> hl.Int32Expression:
        on_diag = (self.get_count(GenotypeType.HOM_REF, GenotypeType.HOM_REF) +
                   self.get_count(GenotypeType.HET, GenotypeType.HET) +
                   self.get_count(GenotypeType.HOM_ALT, GenotypeType.HOM_ALT))
        return on_diag

    def n_discordant(self) -> hl.Int64Expression:
        n_discordant = self.n_total() - self.n_concordant()
        return n_discordant

    def n_total(self) -> hl.NumericExpression:
        return hl.sum(hl.map(lambda x: hl.sum(x[1:]), self._data[2:]))

    def het_to_ref(self) -> hl.Int32Expression:
        return self.get_count(GenotypeType.HET, GenotypeType.HOM_REF)

    def het_to_het(self) -> hl.Int32Expression:
        return self.get_count(GenotypeType.HET, GenotypeType.HET)

    def het_to_alt(self) -> hl.Int32Expression:
        return self.get_count(GenotypeType.HET, GenotypeType.HOM_ALT)

    def n_hets(self) -> hl.NumericExpression:
        return hl.sum(self._data[GenotypeType.HET.value])

    def n_hom_alts(self) -> hl.NumericExpression:
        return hl.sum(self._data[GenotypeType.HOM_ALT.value])

    def flattened_concordance_array(self) -> Dict[str, hl.Expression]:
        return {
            'MISSING_MISSING': self.get_count(GenotypeType.MISSING, GenotypeType.MISSING),
            'MISSING_NOCALL': self.get_count(GenotypeType.MISSING, GenotypeType.NO_CALL),
            'MISSING_HOMREF': self.get_count(GenotypeType.MISSING, GenotypeType.HOM_REF),
            'MISSING_HET': self.get_count(GenotypeType.MISSING, GenotypeType.HET),
            'MISSING_HOMALT': self.get_count(GenotypeType.MISSING, GenotypeType.HOM_ALT),
            'NOCALL_MISSING': self.get_count(GenotypeType.NO_CALL, GenotypeType.MISSING),
            'NOCALL_NOCALL': self.get_count(GenotypeType.NO_CALL, GenotypeType.NO_CALL),
            'NOCALL_HOMREF': self.get_count(GenotypeType.NO_CALL, GenotypeType.HOM_REF),
            'NOCALL_HET': self.get_count(GenotypeType.NO_CALL, GenotypeType.HET),
            'NOCALL_HOMALT': self.get_count(GenotypeType.NO_CALL, GenotypeType.HOM_ALT),
            'HOMREF_MISSING': self.get_count(GenotypeType.HOM_REF, GenotypeType.MISSING),
            'HOMREF_NOCALL': self.get_count(GenotypeType.HOM_REF, GenotypeType.NO_CALL),
            'HOMREF_HOMREF': self.get_count(GenotypeType.HOM_REF, GenotypeType.HOM_REF),
            'HOMREF_HET': self.get_count(GenotypeType.HOM_REF, GenotypeType.HET),
            'HOMREF_HOMALT': self.get_count(GenotypeType.HOM_REF, GenotypeType.HOM_ALT),
            'HET_MISSING': self.get_count(GenotypeType.HET, GenotypeType.MISSING),
            'HET_NOCALL': self.get_count(GenotypeType.HET, GenotypeType.NO_CALL),
            'HET_HOMREF': self.get_count(GenotypeType.HET, GenotypeType.HOM_REF),
            'HET_HET': self.get_count(GenotypeType.HET, GenotypeType.HET),
            'HET_HOMALT': self.get_count(GenotypeType.HET, GenotypeType.HOM_ALT),
            'HOMALT_MISSING': self.get_count(GenotypeType.HOM_ALT, GenotypeType.MISSING),
            'HOMALT_NOCALL': self.get_count(GenotypeType.HOM_ALT, GenotypeType.NO_CALL),
            'HOMALT_HOMREF': self.get_count(GenotypeType.HOM_ALT, GenotypeType.HOM_REF),
            'HOMALT_HET': self.get_count(GenotypeType.HOM_ALT, GenotypeType.HET),
            'HOMALT_HOMALT': self.get_count(GenotypeType.HOM_ALT, GenotypeType.HOM_ALT),
        }

    def result(self) -> Dict[str, hl.Expression]:
        return {
            'concordance_rate': hl.or_missing(hl.is_defined(self.concordance_rate()), hl.format('%.2f', self.concordance_rate() * 100)),
            'nonref_concordance_rate': hl.or_missing(hl.is_defined(self.nonref_concordance()), hl.format('%.2f', self.nonref_concordance() * 100)),
            'f1_score': hl.or_missing(hl.is_defined(self.f1_score()), hl.format('%.2f', self.f1_score() * 100)),
            'n_total': self.n_total(),
            'n_concordant': self.n_concordant(),
            'n_discordant': self.n_discordant(),
            'n_hets': self.n_hets(),
            'n_hom_alts': self.n_hom_alts(),
            'het_to_hom_alt': self.get_count(GenotypeType.HET, GenotypeType.HOM_ALT),
            'hom_ref_to_het': self.get_count(GenotypeType.HOM_REF, GenotypeType.HET),
            'het_to_hom_ref': self.get_count(GenotypeType.HET, GenotypeType.HOM_REF),
        }


class ConcordanceView:
    #: A list of all possible fields in the view.
    fields = [
        'concordance_rate',
        'nonref_concordance_rate',
        'f1_score',
        'n_total',
        'n_concordant',
        'n_discordant,'
        'n_hets',
        'n_hom_alts',
        'het_to_hom_alt',
        'hom_ref_to_het',
        'het_to_hom_ref',
        'MISSING_MISSING',
        'MISSING_NOCALL',
        'MISSING_HOMREF',
        'MISSING_HET',
        'MISSING_HOMALT',
        'NOCALL_MISSING',
        'NOCALL_NOCALL',
        'NOCALL_HOMREF',
        'NOCALL_HET',
        'NOCALL_HOMALT',
        'HOMREF_MISSING',
        'HOMREF_NOCALL',
        'HOMREF_HOMREF',
        'HOMREF_HET',
        'HOMREF_HOMALT',
        'HET_MISSING',
        'HET_NOCALL',
        'HET_HOMREF',
        'HET_HET',
        'HET_HOMALT',
        'HOMALT_MISSING',
        'HOMALT_NOCALL',
        'HOMALT_HOMREF',
        'HOMALT_HET',
        'HOMALT_HOMALT',
    ]

    def __init__(self,
                 *,
                 table: hl.Table,
                 group_vars: Dict[str, hl.HailType],
                 agg_vars: Dict[str, hl.HailType],
                 ordering: Dict[str, List[str]]):
        self._table = table
        self._group_vars = group_vars
        self._agg_vars = agg_vars
        self._ordering = ordering

        self._df: Optional[pd.DataFrame] = None
        self._key_fields: Optional[List[str]] = None

    def _to_df(self) -> Tuple[pd.DataFrame, List[str]]:
        if self._df is not None:
            return (self._df, self._key_fields)

        key_fields = list(self._table.key)

        df = self._table.to_pandas()
        df = df.dropna(subset=['group_var', 'group_var_name', 'agg_var', 'agg_var_name']).copy()

        def get_ordering(var_name):
            key = tuple(var_name) if isinstance(var_name, list) else var_name
            return self._ordering.get(key, [])

        def cast_value(value, cast_type):
            try:
                return cast_type(value)
            except (ValueError, TypeError):
                return value

        def cast_variable(row):
            group_name = row['group_var_name']
            group_type = get_cast_type(self._group_vars[group_name])
            if pd.notna(row['group_var']):
                row['group_var'] = cast_value(row['group_var'], group_type)

            agg_name = row['agg_var_name']
            agg_type = get_cast_type(self._agg_vars[agg_name])
            if pd.notna(row['agg_var']):
                row['agg_var'] = cast_value(row['agg_var'], agg_type)

            return row

        dfs = []
        for group in self._group_vars.keys():
            for agg in self._agg_vars.keys():
                tmp = df[(df['group_var_name'] == group) & (df['agg_var_name'] == agg)].copy()

                tmp = tmp.apply(cast_variable, axis=1)

                if group in self._ordering:
                    tmp['group_var'] = pd.Categorical(tmp['group_var'], categories=get_ordering(group), ordered=True)

                if agg in self._ordering:
                    tmp['agg_var'] = pd.Categorical(tmp['agg_var'], categories=get_ordering(agg), ordered=True)

                dfs.append(tmp)

        df = pd.concat(dfs, ignore_index=True)

        df = df.sort_values(by=['group_var_name', 'group_var', 'agg_var_name', 'agg_var'])

        self._df = df
        self._key_fields = key_fields

        return (df, key_fields)

    def result(self) -> pd.DataFrame:
        """Generate a pandas DataFrame of the underlying view.

        Examples:
            >>> view = conc_table.get_view(group_names=['exo_mac'], agg_names=['DP', 'GQ'], ordering={'exo_mac': ['1', '2-5', '6-10', '10+']})
            >>> df = view.result()
            >>> df.head(5)
        """
        return self._to_df()[0]

    def _split_results(self, fields: Optional[Dict[str, str]] = None) -> Dict[Tuple[str, str], pd.DataFrame]:
        results = {}
        for group in self._group_vars.keys():
            for agg in self._agg_vars.keys():
                df, key_fields = self._to_df()

                df = df[(df['group_var_name'] == group) & (df['agg_var_name'] == agg)].copy()

                df[group] = df['group_var']
                df[agg] = df['agg_var']
                df = df.drop(columns=['group_var_name', 'group_var', 'agg_var_name', 'agg_var', 'concordance'])

                df = df[key_fields + [group, agg] + list(fields.keys())].copy()

                if group in self._ordering:
                    df[group] = pd.Categorical(df[group], categories=self._ordering[group], ordered=True)

                if agg in self._ordering:
                    df[agg] = pd.Categorical(df[agg], categories=self._ordering[agg], ordered=True)

                if fields is not None:
                    df.rename(columns=fields, inplace=True)

                results[(group, agg)] = df.copy()
        return results

    def export(self, path: str, fields: Optional[Dict[str, str]] = None, **kwargs):
        """Export the underlying Pandas dataframe to a file.

        Notes:
            The `**kwargs` are passed to the "to_csv" method of a pandas DataFrame.

        Examples:

            >>> view = conc_table.get_view(group_names=['exo_mac'], agg_names=['DP', 'GQ'], ordering={'exo_mac': ['1', '2-5', '6-10', '10+']})
            >>> view.export('gs://my-bucket/results.tsv', sep="\t")

            The contents of `results.tsv` contains all combinations of group variables and aggregation variables.

        Args:
            path (str): Path to export results to.
            fields (Optional[Dict[str, str]]): Fields to export along with the new field name.
            kwargs (dict): Optional kwargs to pass to "pd.DataFrame.to_csv()"
        """
        df, key_fields = self._to_df()

        if fields is not None:
            for new_name, old_name in fields.items():
                df[new_name] = df[old_name]
            df = df[key_fields + list(fields.keys())]

        df = df.drop(columns=['concordance'])

        with NamedTemporaryFile(delete=False, mode='w', newline='') as temp_file:
            df.to_csv(temp_file, **kwargs)
            temp_path = temp_file.name

        hfs.copy(temp_path, path)

    def export_all(self, output_dir: str, *args, fields: Optional[Dict[str, str]] = None, **kwargs):
        """Export all possible combinations of group and agg variables.

        Examples:

            >>> view = conc_table.get_view(group_names=['exo_mac'], agg_names=['DP', 'GQ'], ordering={'exo_mac': ['1', '2-5', '6-10', '10+']})
            >>> view.export('gs://my-bucket/results/', sep="\t")

            There is one table exported for each combination of grouping variable and aggregation variable. The files have the structure "GROUP_AGG.tsv".
            For example, the view above would generate the following files:

            - exo_mac_DP.tsv
            - exo_mac_GQ.tsv

        Args:
            output_dir (str): Directory to write all files to.
            fields (Optional[Dict[str, str]]): Fields to export along with the new field name
        """
        output_dir = output_dir.rstrip('/') + '/'
        kwargs['sep'] = '\t'

        dfs = self._split_results(fields=fields)

        for (group, agg), df in dfs.items():
            path = f'{output_dir}{group}_{agg}.tsv'
            with NamedTemporaryFile(delete=False, mode='w', newline='') as temp_file:
                df.to_csv(temp_file, *args, **kwargs)
                temp_path = temp_file.name

            hfs.copy(temp_path, path)

    def _plot_base(self, df: pd.DataFrame, *, statistic: Statistic = Statistic.NONREF_CONCORDANCE) -> ggplot:
        df = df[df[statistic.value].notnull()].copy()
        df[statistic.value] = df[statistic.value].astype(float)

        group = df['group_var_name'].tolist()[0]
        n_aggs = len(self._agg_vars)

        plot = (
                ggplot(df, aes(x="agg_var", y=statistic.value, color="group_var"))
                + geom_point(aes(size="n_total"), alpha=0.6, show_legend=False)
                + geom_line(aes(group="group_var"), size=1.5)
                + scale_size_continuous(name="Number of Genotypes")
                + theme_minimal()
                + theme(
            figure_size=(3 * n_aggs, 3),
            legend_position="right",
            axis_text_x=element_text(rotation=90, hjust=1, size=12),
            plot_title=element_text(size=18),
            axis_title_x=element_text(size=14),
            axis_title_y=element_text(size=14),
            axis_text_y=element_text(size=12),
            legend_title=element_text(size=14),
            legend_text=element_text(size=12),
            panel_background=element_rect(fill='white'),
            plot_background=element_rect(fill='white')
        )
                + coord_fixed(ratio=1.5)
                + labs(title="Concordance", y=statistic.value, x="")
                + facet_grid("~ agg_var_name", scales="free_x")
                + guides(color=guide_legend(title=group))
        )

        group_order = self._ordering.get(group)
        if group_order is not None:
            plot = plot + scale_color_discrete(limits=group_order)

        return plot

    def plot(self, *, path: Optional[str] = None, statistic: Statistic = Statistic.NONREF_CONCORDANCE, **kwargs) -> ggplot:
        """Create a plot of the data that features a facet per aggregation variable in the view.

        Notes:
            Only implemented for one grouping variable in the view. Otherwise, use `plot_all` or create a new view with one grouping variable.

        Examples:
            >>> view = conc_table.get_view(group_names=['exo_mac'], agg_names=['DP', 'GQ'], ordering={'exo_mac': ['1', '2-5', '6-10', '10+']})
            >>> p = view.plot()
            >>> p.save(filename='plot.png')

        Args:
              path (Optional[str]): An optional path to write the figure to.
              statistic (Statistic): A statistic to use as the plotting variable. One of Statistic.NONREF_CONCORDANCE and Statistic.F1_SCORE.
              kwargs: Pass-through arguments to plotnine.save

        Returns
            plotnine.ggplot: A plot object
        """
        if len(self._group_vars) == 1:
            df = self._to_df()[0]
            p = self._plot_base(df, statistic=statistic)
            if path is not None:
                with NamedTemporaryFile(delete=False, suffix='.png') as temp_file:
                    temp_file_name = temp_file.name

                p.draw()
                p.save(temp_file_name, transparent=False, **kwargs)

                hfs.copy(temp_file_name, path)
                os.remove(temp_file_name)
            return p
        raise NotImplementedError('no plotting support for multiple group variables')

    def plot_all(self,
                 *,
                 out_dir: Optional[str] = None,
                 statistic: Statistic = Statistic.NONREF_CONCORDANCE,
                 **kwargs) -> Dict[str, ggplot]:
        """Create plots that feature a facet per aggregation variable in the view.

        Args:
              out_dir (Optional[str]): An optional directory to save plots to
              statistic (Optional[Statistic]): A statistic to use as the plotting variable.
              kwargs: Keyword arguments to pass through to plotnine.save

        Returns
            Dict[str, plotnine.ggplot]: A dictionary mapping a group name to the corresponding plot object.
        """
        plots = {}
        for group in self._group_vars.keys():
            df = self._to_df()[0]
            df = df[(df['group_var_name'] == group)].copy()

            p = self._plot_base(df, statistic=statistic)

            if out_dir is not None:
                out_dir = out_dir.rstrip('/') + '/'
                with NamedTemporaryFile(delete=False, suffix='.png') as temp_file:
                    temp_file_name = temp_file.name

                p.draw()
                p.save(temp_file_name, transparent=False, **kwargs)

                hfs.copy(temp_file_name, f'{out_dir}{group}.png')
                os.remove(temp_file_name)

            plots[group] = p

        return plots


class ConcordanceTable:
    """Result of running a concordance operation on a JointCallset."""

    @staticmethod
    def load(path: str):
        """Path to a Hail Table containing the raw concordance output.

        Examples:
            >>> conc_table = ConcordanceTable.load('gs://my-bucket/global-conc.ht')

        Args:
            path (str): The path of the Hail Table to import.

        Returns:
            ConcordanceTable: A class for interacting with the results of a concordance operation.
        """
        t = hl.read_table(path)
        return ConcordanceTable(t)

    @staticmethod
    def _from_concordance_output(t: hl.Table) -> 'ConcordanceTable':
        t = t.annotate(concordance=hl.array(t.concordance))
        t = t.explode('concordance')
        t = t.annotate(group_var_name=t.concordance[0], concordance=hl.array(t.concordance[1]))
        t = t.explode('concordance')
        t = t.annotate(group_var=t.concordance[0], concordance=hl.array(t.concordance[1]))
        t = t.explode('concordance')
        t = t.annotate(agg_var_name=t.concordance[0], concordance=hl.array(t.concordance[1]))
        t = t.explode('concordance')
        t = t.annotate(agg_var=t.concordance[0], concordance=t.concordance[1])
        t = t.annotate(n_variants=t.concordance.n_variants[0], n_samples=t.concordance.n_samples, concordance=t.concordance.concordance)
        t = t.annotate(**ConcordanceArray(t.concordance).result(),
                       **ConcordanceArray(t.concordance).flattened_concordance_array())
        return ConcordanceTable(t)

    def __init__(self, table: hl.Table):
        self._table = table
        self._group_types = None
        self._agg_types = None
        self._views = []

    def close(self):
        """Unpersist any cached views."""
        self._table.unpersist()
        for view in self._views:
            view.unpersist()

    def write(self, output_path: str, overwrite: bool = True):
        """Write the underlying data to a Hail Table.

        Example:
            >>> concordance_table.write('/tmp/my-output/my-data.ht')

        Args:
            output_path (str): Output path to write the underlying data to.
            overwrite (bool): Overwrite any existing files.
        """
        self._table.write(output_path, overwrite=overwrite)
        self._table = hl.read_table(output_path)

    def group_vars(self) -> Dict[str, hl.HailType]:
        """A mapping of the group by variable names and their types.

        Example:
            >>> concordance_table.group_vars()
            {'exome_mac': hl.tint32, 'exome_maf': hl.tfloat64}

        Returns:
            Dict[str, HailType]: A dictionary with the group by variable names and the corresponding types of their values.
        """
        if self._group_types is None:
            group_types = self._table.group_types.collect()[0]
            self._group_types = {name: hl.dtype(typ) for name, typ in group_types.items()}
        return self._group_types

    def agg_vars(self) -> Dict[str, hl.HailType]:
        """A mapping of the aggregation variable names and their types.

        Example:
            >>> concordance_table.group_vars()
            {'DP': hl.tint32, 'INFO': hl.tfloat64}

        Returns:
            Dict[str, HailType]: A dictionary with the aggregation variable names and the corresponding types of their values.
        """
        if self._agg_types is None:
            agg_types = self._table.agg_types.collect()[0]
            self._agg_types = {name: hl.dtype(typ) for name, typ in agg_types.items()}
        return self._agg_types

    def get_view(self,
                 *,
                 key: Optional[Dict[str, Any]] = None,
                 group_names: Optional[List[str]] = None,
                 agg_names: Optional[List[str]] = None,
                 ordering: Optional[Dict[str, List[str]]] = None) -> ConcordanceView:
        """Subset the concordance results for a particular key of the concordance results.

       Examples:

           >>> global_view = glob_conc_table.get_view(group_names=['exo_mac'], agg_names=['DP', 'GQ'], ordering={'exo_mac': ['1', '2-5', '6-10', '10+']})
           >>> global_view.plot()
           >>> global_view.export('gs://my-bucket/my-file.tsv', sep='\t')

           >>> variant_view = variant_conc_table.get_view(key=dict(locus=hl.Locus('chr1', 3423242, reference_genome='GRCh38'), alleles=['A', 'G']))

           >>> sample_view = sample_conc_table.get_view(key=dict(s="NA12878"))

       Args:
           key (Optional[Dict[str, Any]]): The key to subset the data to. For a global concordance result, this is None. For a variant concordance result, this is the locus and alleles. For sample concordance, this is the sample ID.
           group_names (Optional[List[str]]): The group by variables to select. The default is all variants and all samples or no grouping.
           agg_names (Optional[List[str]]): The aggregation variables to select. The default is all genotypes.
           ordering (Optional[Dict[str, List[str]]]): A mapping from group name or agg name to an ordered list of the factors for that variable.

       Returns:
           ConcordanceView: A slice of the data that can be exported or plotted.
        """
        t = self._table

        key = key or {}
        group_names = group_names or list(self.group_vars().keys())
        agg_names = agg_names or list(self.agg_vars().keys())

        filter_expr = hl.bool(True)

        if key is not None:
            for key_name, value in key.items():
                filter_expr &= (t[key_name] == value)

        if group_names is not None:
            filter_expr &= hl.literal(group_names).contains(t['group_var_name'])

        if agg_names is not None:
            filter_expr &= hl.literal(agg_names).contains(t['agg_var_name'])

        t = t.filter(filter_expr)

        group_vars = {k: v for k, v in self.group_vars().items() if k in group_names}
        agg_vars = {k: v for k, v in self.agg_vars().items() if k in agg_names}

        t = t.persist()
        self._views.append(t)

        return ConcordanceView(table=t, group_vars=group_vars, agg_vars=agg_vars, ordering=ordering or {})


def concordance_array(joined):
    def get_idx(struct):
        return hl.if_else(hl.is_missing(struct), 0, hl.coalesce(2 + struct.GT.n_alt_alleles(), 1))

    counter = hl.agg.group_by(get_idx(joined.left_entry) + 5 * get_idx(joined.right_entry), hl.agg.count())

    default_value = hl.int64(0)

    return hl.range(0, 5).map(lambda i: hl.range(0, 5).map(lambda j: counter.get(i + 5 * j, default_value)))


def variant_count(joined):
    reference_genome = joined.locus.dtype.reference_genome
    first_variant = hl.locus(reference_genome.contigs[0], 1, reference_genome=reference_genome)

    return hl.agg.fold(hl.tuple([0, first_variant]),
                       lambda accum: hl.if_else(joined.left_row.locus != accum[1], hl.tuple([accum[0] + 1, joined.left_row.locus]), hl.tuple([accum[0], joined.left_row.locus])),
                       lambda comb_left, comb_right: hl.tuple([comb_left[0] + comb_right[0], comb_right[1]]))


def sample_count(joined):
    return hl.len(hl.agg.collect_as_set(joined.col_key))


def global_concordance(exo: hl.MatrixTable,
                       imp: hl.MatrixTable,
                       *,
                       exo_row_group_by_var: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                       imp_row_group_by_var: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                       exo_col_group_by_var: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                       imp_col_group_by_var: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                       exo_row_agg_fields: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                       imp_row_agg_fields: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                       exo_col_agg_fields: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                       imp_col_agg_fields: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                       exo_entry_agg_fields: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                       imp_entry_agg_fields: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                       join_type: JoinType = JoinType.inner) -> ConcordanceTable:
    joined = full_outer_join(exo, imp, join_type)

    group_names = [ALL_GROUP,
                   *list(exo_row_group_by_var.keys()),
                   *list(imp_row_group_by_var.keys()),
                   *list(exo_col_group_by_var.keys()),
                   *list(imp_col_group_by_var.keys())]

    group_info = [*[(ALL_GROUP, hl.tstr)],
                  *[(f(joined.left_row), f(joined.left_row).dtype) for name, f in exo_row_group_by_var.items()],
                  *[(f(joined.right_row), f(joined.right_row).dtype) for name, f in imp_row_group_by_var.items()],
                  *[(f(joined.left_row), f(joined.left_row).dtype) for name, f in exo_col_group_by_var.items()],
                  *[(f(joined.right_row), f(joined.right_row).dtype) for name, f in imp_col_group_by_var.items()]]

    group_values = hl.array([hl.str(group) for group, _ in group_info])
    group_types = {group_name: str(typ) for group_name, (_, typ) in zip(group_names, group_info)}

    agg_names = [ALL_AGG,
                 *exo_row_agg_fields,
                 *imp_row_agg_fields,
                 *exo_col_agg_fields,
                 *imp_col_agg_fields,
                 *exo_entry_agg_fields,
                 *imp_entry_agg_fields]

    agg_info = [*[(ALL_AGG, hl.tstr)],
                *[(f(joined.left_row), f(joined.left_row).dtype) for f in exo_row_agg_fields.values()],
                *[(f(joined.right_row), f(joined.right_row).dtype) for f in imp_row_agg_fields.values()],
                *[(f(joined.left_col), f(joined.left_col).dtype) for f in exo_col_agg_fields.values()],
                *[(f(joined.right_col), f(joined.right_col).dtype) for f in imp_col_agg_fields.values()],
                *[(f(joined.left_entry), f(joined.left_entry).dtype) for f in exo_entry_agg_fields.values()],
                *[(f(joined.right_entry), f(joined.right_entry).dtype)  for f in imp_entry_agg_fields.values()]]

    agg_values = hl.array([hl.str(agg) for agg, _ in agg_info])
    agg_types = {agg_name: str(typ) for agg_name, (_, typ) in zip(agg_names, agg_info)}

    result_agg = hl.struct(n_samples=sample_count(joined),
                           n_variants=variant_count(joined),
                           concordance=concordance_array(joined))

    aggr = hl.agg.array_agg(lambda key: hl.agg.group_by(key, result_agg), agg_values)

    result = joined.aggregate_entries(hl.struct(
        concordance=hl.dict(hl.zip(hl.literal(group_names), hl.agg.array_agg(lambda key: hl.agg.group_by(key, aggr), group_values))).map_values(
            lambda grouped_d: grouped_d.map_values(lambda arr: hl.dict(hl.zip(hl.literal(agg_names), arr)))
        )
    ), _localize=False)

    t = hl.utils.range_table(1)
    t = t.annotate(**result)
    t = t.annotate_globals(group_types=group_types, agg_types=agg_types)
    t = t.key_by()
    t = t.drop('idx')

    return ConcordanceTable._from_concordance_output(t)


def variant_concordance(exo: hl.MatrixTable,
                        imp: hl.MatrixTable,
                        *,
                        exo_col_group_by_var: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                        imp_col_group_by_var: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                        exo_col_agg_fields: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                        imp_col_agg_fields: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                        exo_entry_agg_fields: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                        imp_entry_agg_fields: Dict[str, Callable[[hl.StructExpression], hl.Expression]],
                        join_type: JoinType = JoinType.inner) -> ConcordanceTable:
    joined = full_outer_join(exo, imp, join_type)

    group_names = [ALL_GROUP,
                   *list(exo_col_group_by_var.keys()),
                   *list(imp_col_group_by_var.keys())]

    group_info = [*[(ALL_GROUP, hl.tstr)],
                  *[(f(joined.left_row), f(joined.left_row).dtype) for name, f in exo_col_group_by_var.items()],
                  *[(f(joined.right_row), f(joined.right_row).dtype) for name, f in imp_col_group_by_var.items()]]

    group_values = hl.array([hl.str(group) for group, _ in group_info])
    group_types = {group_name: str(typ) for group_name, (_, typ) in zip(group_names, group_info)}

    agg_names = [ALL_AGG,
                 *exo_col_agg_fields,
                 *imp_col_agg_fields,
                 *exo_entry_agg_fields,
                 *imp_entry_agg_fields]

    agg_info = [*[(ALL_AGG, hl.tstr)],
                *[(f(joined.left_col), f(joined.left_col).dtype) for f in exo_col_agg_fields.values()],
                *[(f(joined.right_col), f(joined.right_col).dtype) for f in imp_col_agg_fields.values()],
                *[(f(joined.left_entry), f(joined.left_entry).dtype) for f in exo_entry_agg_fields.values()],
                *[(f(joined.right_entry), f(joined.right_entry).dtype) for f in imp_entry_agg_fields.values()]]

    agg_values = hl.array([hl.str(agg) for agg, _ in agg_info])
    agg_types = {agg_name: str(typ) for agg_name, (_, typ) in zip(agg_names, agg_info)}

    result_agg = hl.struct(n_samples=sample_count(joined),
                           n_variants=variant_count(joined),
                           concordance=concordance_array(joined))

    aggr = hl.agg.array_agg(lambda key: hl.agg.group_by(key, result_agg), agg_values)

    result = joined.annotate_rows(concordance=hl.dict(hl.zip(hl.literal(group_names), hl.agg.array_agg(lambda key: hl.agg.group_by(key, aggr), group_values))).map_values(
            lambda grouped_d: grouped_d.map_values(lambda arr: hl.dict(hl.zip(hl.literal(agg_names), arr)))))

    result = result.annotate_globals(group_types=group_types, agg_types=agg_types)
    result = result.rows()
    result = result.select('concordance')
    result = result.key_by(contig=result.locus.contig, position=result.locus.position, ref=result.alleles[0], alt=hl.delimit(result.alleles[1:]))
    result = result.drop('locus', 'alleles')
    return ConcordanceTable._from_concordance_output(result)


def sample_concordance(exo: hl.MatrixTable, imp: hl.MatrixTable, *, join_type: JoinType = JoinType.inner) -> ConcordanceTable:
    joined = full_outer_join(exo, imp, join_type)

    result_agg = hl.struct(n_samples=sample_count(joined),
                           n_variants=variant_count(joined),
                           concordance=concordance_array(joined))

    result = joined.annotate_cols(concordance=hl.dict([(ALL_GROUP, hl.dict([(ALL_GROUP, hl.dict([(ALL_AGG, hl.dict([(ALL_AGG, result_agg)]))]))]))]))
    result = result.annotate_globals(group_types={ALL_GROUP: str(hl.tstr)}, agg_types={ALL_AGG: str(hl.tstr)})
    result = result.cols().select('concordance')

    return ConcordanceTable._from_concordance_output(result)
