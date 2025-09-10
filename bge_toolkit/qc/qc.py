from typing import Dict, List, Optional, Union
import logging
import hail as hl
from bge_toolkit.common.utils import setup_logger
from hailtop.aiocloud.aiogoogle import GCSRequesterPaysConfiguration
import hailtop.fs as hfs
from gnomad.sample_qc.filtering import compute_stratified_metrics_filter
from gnomad.sample_qc.ancestry import apply_onnx_classification_model, assign_population_pcs

import pandas as pd

from plotnine import (ggplot, aes, geom_point, scale_color_manual, theme_minimal, labs, theme,
                      element_text, coord_fixed, element_rect, geom_density, geom_jitter,
                      geom_boxplot, coord_flip)

import onnx

log = logging.getLogger()
log.setLevel(logging.INFO)


filtering_qc_metrics = ['r_ti_tv',
                        'n_singleton',
                        'n_insertion',
                        'n_deletion',
                        'n_transition',
                        'n_transversion',
                        'r_het_hom_var',
                        'r_insertion_deletion',
                        ]


def filter_to_high_quality_exome_intervals(dataset: Union[hl.MatrixTable, hl.vds.VariantDataset],
                                           low_complexity_regions: hl.Table,
                                           exome_regions: hl.Table,
                                           passing_variants: hl.Table,
                                           is_gatk: bool = True) -> hl.MatrixTable:
    """Filter a MatrixTable to the exome with no low complexity regions.

    Args:
        dataset (MatrixTable, VariantDataset): A Hail MatrixTable or VariantDataset to filter the sites.
        low_complexity_regions (Table): A Table of intervals to filter out. Use `import_locus_intervals` from a BED file to generate the table.
        exome_regions (Table): A Table of intervals to filter to. Use `import_locus_intervals` from a BED file to generate the table.
        passing_variants (Table): A Table of variants to filter to that pass VQSR or filters
        is_gatk (bool): True if the dataset was called using GATK.
    Result:
        A filtered MatrixTable to only exome regions without the low complexity regions of the genome.
    """
    log = setup_logger()

    if isinstance(dataset, hl.MatrixTable):
        dataset = dataset.filter_rows(hl.is_missing(low_complexity_regions[dataset.locus]))
        dataset = dataset.filter_rows(hl.is_defined(exome_regions[dataset.locus]))
        dataset = dataset.filter_rows(hl.is_defined(passing_variants[dataset.locus]))
    else:
        assert isinstance(dataset, hl.vds.VariantDataset)
        log.info('filtering intervals for exome regions')
        variant_data = dataset.variant_data
        variant_data = variant_data.filter_rows(hl.is_defined(exome_regions[variant_data.locus]))
        dataset = hl.vds.VariantDataset(dataset.reference_data, variant_data)
        # dataset = hl.vds.filter_intervals(dataset, exome_regions, keep=True)

        log.info('filtering intervals for LCR')
        variant_data = dataset.variant_data
        variant_data = variant_data.filter_rows(hl.is_missing(low_complexity_regions[variant_data.locus]))
        dataset = hl.vds.VariantDataset(dataset.reference_data, variant_data)
        # dataset = hl.vds.filter_intervals(dataset, low_complexity_regions, keep=False)

        log.info('filtering variants to passing variants')
        variant_data = dataset.variant_data
        variant_data = variant_data.semi_join_rows(passing_variants)
        dataset = hl.vds.VariantDataset(dataset.reference_data, variant_data)
        # dataset = hl.vds.filter_variants(dataset, passing_variants, keep=True)
    return dataset


def filter_contamination_rate(*,
                              mt: hl.MatrixTable,
                              charr_thresh: float = 0.05,
                              min_af: float = 0.05,
                              max_af: float = 0.95,
                              min_dp: int = 10,
                              max_dp: int = 100,
                              min_gq: int = 20,
                              ref_AF: Optional[hl.Float64Expression] = None) -> hl.Table:
    """Identify samples with high contamination rates.

    Uses the Hail ``compute_charr`` function to compute the estimate of contamination rate.

    The result is a Hail Table with the following fields:
      * mt.col_key: The column key of the input MatrixTable
      * contamination:
        * charr: The CHARR statistic.
        * is_passing: ``charr`` is less than the ``charr_thresh``

    Args:
        mt (MatrixTable): A Hail MatrixTable to compute contamination statistics from.
        charr_thresh (float): The threshold for which a sample is considered to be passing.
        min_af (float): The minimum allele frequency when using the ``compute_charr`` function in Hail.
        max_af (float): The maximum allele frequency when using the ``compute_charr`` function in Hail.
        min_dp (int): The minimum depth when using the ``compute_charr`` function in Hail.
        max_dp (int): The maximum depth when using the ``compute_charr`` function in Hail.
        min_gq (int): The minimum GQ when using the ``compute_charr`` function in Hail.
        ref_AF (Float64Expression): A float row field on the MatrixTable with the reference allele frequency when using the ``compute_charr`` function in Hail.

    Returns:
        A Hail Table with the contamination statistic in a new column "contamination" that is a struct with one field "charr".
    """
    charr_t = hl.compute_charr(mt, min_af=min_af, max_af=max_af, min_dp=min_dp, max_dp=max_dp, min_gq=min_gq, ref_AF=ref_AF)
    t = charr_t.select(contamination=hl.struct(charr=charr_t.charr, is_passing=charr_t.charr <= charr_thresh))
    return t


def filter_chimeric_reads(*,
                          chimera_rate: hl.NumericExpression,
                          threshold: float) -> hl.Table:
    """Determine whether a sample passes based on chimeric read percentages.

    The result is a Hail Table with the following fields:
      * ds.col_key: The column key of the input expression.
      * chimera_reads:
        * chimera_rate: The rate of chimera reads.
        * is_passing: The rate of chimera reads is below ``threshold``.

    Args:
        chimera_rate (NumericExpression): A numeric expression containing the chimeric read rate.
        threshold (float): The maximum threshold for which to consider a sample passing.

    Returns:
        A Hail Table with the chimeric read data along with an annotation for whether it is passing.
    """
    name, chimera_rate = chimera_rate._to_relational_preserving_rows_and_cols('chimera_rate')
    chimera_rate = chimera_rate.select(chimera_reads=hl.struct(chimera_rate=chimera_rate[name],
                                                               is_passing=chimera_rate[name] < threshold))
    return chimera_rate


def filter_genotypes(mt: hl.MatrixTable,
                     is_gatk: bool,
                     *,
                     dp_thresh: int = 10,
                     gq_thresh: int = 20,
                     ab_thresh: float = 0.2
                     ) -> hl.MatrixTable:
    """Filter genotypes that are of poor quality from an exome dataset.

    If the dataset was called using DRAGEN, ``DP`` is the sum of ``AD`` and only heterozygote
    and homozygous alternate calls are subject to the filter on read depth.

    Args:
        mt (MatrixTable): A MatrixTable to filter the genotypes of.
        is_gatk (bool): True if the dataset was called using GATK.
        dp_thresh (int): The minimum DP for a kept genotype. DP is computed as the sum of the "AD" field.
        gq_thresh (int): The minimum GQ threshold for a kept genotype.
        ab_thresh (float): The allowable allelic balance range for a kept heterozygote genotype.

    Returns:
        A MatrixTable with poor quality genotypes removed.
    """
    mt = mt.filter_entries(hl.is_defined(mt.GT))

    if is_gatk:
        mt = mt.filter_entries(mt.DP >= dp_thresh)
    else:
        # DRAGEN data
        mt = mt.annotate_entries(DP=hl.sum(mt.AD))

        dp_filter = (hl.case(missing_false=True)
                     .when(mt.GT.n_alt_alleles() > 0, mt.DP >= dp_thresh)
                     .when(mt.GT.n_alt_alleles() == 0, True)
                     .or_missing())

        mt = mt.filter_entries(dp_filter)

    mt = mt.filter_entries(mt.GQ >= gq_thresh)

    alt_ab = mt.AD[1] / mt.DP
    is_het_ab_within_thresh = (alt_ab >= ab_thresh) & (alt_ab <= (1 - ab_thresh))

    mt = mt.filter_entries(~mt.GT.is_het() | (mt.GT.is_het() & is_het_ab_within_thresh))

    return mt


def filter_to_high_quality_common_sites(mt: hl.MatrixTable,
                                        *,
                                        maf_thresh: float = 0.01,
                                        mac_thresh: int = 10,
                                        call_rate_thresh: float = 0.95) -> hl.MatrixTable:
    """Filter a MatrixTable to common, LD pruned sites.

    Args:
        mt (MatrixTable): The MatrixTable to subset to filter to common, LD pruned sites.
        maf_thresh (float): The minimum minor allele frequency as computed by `hl.variant_qc` in order for a variant to be considered to be common.
        mac_thresh (int): The minimum number of minor alleles observed as computed by `hl.variant_qc` in order for a variant to be considered common.
        call_rate_thresh (float): The minimum call rate of a variant as computed by `hl.variant_qc` in order for a variant to be considered.

    Returns:
        A Hail MatrixTable subsetted to common, LD-pruned sites.
    """
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(hl.min(mt.variant_qc.AF) > maf_thresh)
    mt = mt.filter_rows(mt.variant_qc.call_rate > call_rate_thresh)
    mt = mt.filter_rows(hl.min(mt.variant_qc.AC) >= mac_thresh)

    return mt


def compute_pcs(mt: hl.MatrixTable, *, k: int = 10) -> hl.Table:
    """Compute principal components.

    Only computed from autosomal variants.

    The result is a Hail Table with the following fields:
      * mt.col_key: The column key of the input MatrixTable.
      * pcs:
        * scores: An array with the top ``k`` principal components.

    Args:
        mt (MatrixTable): A Hail MatrixTable to compute principal components from. Should be LD pruned and filtered to high quality common variants.
        k (int): The number of principal components to compute.

    Returns:
        A Hail Table with the principal components.
    """
    mt = mt.filter_rows(mt.locus.in_autosome())
    _, scores, _ = hl.hwe_normalized_pca(mt.GT, k=k)
    scores = scores.select(pcs=hl.struct(scores=scores.scores))
    return scores


def infer_relatedness(mt: hl.MatrixTable,
                      min_individual_maf: float = 0.01,
                      k: Optional[int] = None,
                      scores_expr: Optional[hl.ArrayNumericExpression] = None,
                      min_kinship: Optional[float] = 0.125) -> hl.Table:
    """Identify unrelated samples.

    Only computed for autosomal variants.

    The result is a Hail Table with the following fields:
      * mt.col_key: The column key of the input MatrixTable
      * relatedness:
        * is_passing: An unrelated sample for downstream analysis.

    Args:
        mt (MatrixTable): A MatrixTable to infer relatedness from. This should be a MatrixTable subsetted to high quality, LD pruned common variants.
        min_individual_maf (float): Parameter to the Hail ``pc_relate`` function. The minimum individual-specific minor allele frequency.
          If either individual-specific minor allele frequency for a pair of
          individuals is below this threshold, then the variant will not
          be used to estimate relatedness for the pair.
        k (Optional[int]): Parameter to the Hail ``pc_relate`` function. If set, `k` principal component scores are computed and used.
          Exactly one of `k` and `scores_expr` must be specified. Parameter to the hail ``pc_relate`` function.
        scores_expr (Optional[hl.ArrayNumericExpression]): Parameter to the Hail ``pc_relate`` function. Column-indexed expression of principal component scores, with the same
          source as `call_expr`. All array values must have the same positive length,
          corresponding to the number of principal components, and all scores must
          be non-missing. Exactly one of `k` and `scores_expr` must be specified.
        min_kinship (Optional[float]): Parameter to the Hail ``pc_relate`` function. If set, pairs of samples with kinship lower than `min_kinship` are excluded
          from the results.

    Returns:
        A Hail Table with a flag for the set of unrelated individuals.
    """
    if k is None and scores_expr is None:
        raise ValueError('must specify one of "k" or "scores_expr"')
    if k is not None and scores_expr is not None:
        raise ValueError('must specify only one of "k" or "scores_expr')

    mt = mt.filter_rows(mt.locus.in_autosome())

    pairs = hl.pc_relate(mt.GT,
                         min_individual_maf=min_individual_maf,
                         k=k,
                         scores_expr=scores_expr,
                         statistics='kin',
                         min_kinship=min_kinship)

    unrelated = hl.maximal_independent_set(pairs.i, pairs.j, keep=True)
    unrelated_samples = hl.set(unrelated.node.collect(_localize=False))

    samples = mt.cols().select()
    samples = samples.annotate(relatedness=hl.struct(is_passing=unrelated_samples.contains(samples.key)))
    return samples


def estimate_call_rate_through_gq_stats(mt: hl.MatrixTable,
                                        gq_thresh: int = 20,
                                        fraction_thresh: float = 0.9) -> hl.Table:
    """Flag samples with the proportion of calls with GQ less than a threshold being too low.

    Only computed for autosomal variants.

    The result is a Hail Table with the following fields:
      * mt.col_key: The column key of the input MatrixTable
      * gq_fraction:
        * gq_fraction: The percentage of genotypes with GQ >= ``gq_thresh``.
        * is_passing: ``gq_fraction`` is greater than the ``fraction_thresh``

    Args:
        mt (MatrixTable): The MatrixTable to compute GQ stats from.
        gq_thresh (int): The minimum GQ for a genotype to be considered passing.
        fraction_thresh (float): The minimum rate at which GQ is above the minimum threshold.

    Returns:
        A Hail Table with the struct ``gq_fraction`` and ``include`` in the column ``gq_fraction``.
    """
    mt = mt.filter_rows(mt.locus.in_autosome())

    gq_fraction = hl.agg.fraction(mt.GQ >= gq_thresh)
    mt = mt.annotate_cols(gq_fraction=hl.struct(gq_fraction=gq_fraction, is_passing=gq_fraction >= fraction_thresh))
    return mt.cols().select('gq_fraction')


def impute_and_check_sex(mt: hl.MatrixTable,
                         maf_threshold: float = 0.0,
                         male_fhet_thresh: float = 0.8,
                         female_fhet_thresh: float = 0.2,
                         maf: Optional[hl.NumericExpression] = None,
                         reported_is_female: Optional[hl.BooleanExpression] = None) -> hl.Table:
    """Impute and Check Sex.

    The result is a Hail Table with the following fields:
      * mt.col_key: The column key of the input MatrixTable
      * sex_info:
        * is_female: The imputed sex. ``True`` is for females, ``False`` is for males.
        * is_passing: Either ``reported sex == imputed sex`` or ``True``.
        * sex_check: ``reported sex == imputed sex`` or ``Null``.
        * f_stat: Inbreeding coefficient on the non-PAR X chromosome.
        * n_called: Number of genotypes considered.
        * expected_homs: Expected number of homozygotes.
        * observed_homs: Observed number of heterozygotes.

    Args:
        mt (MatrixTable): A Hail MatrixTable to impute sex from.
        maf_threshold (float): The minor allele frequency threshold in order for a variant to be included in the analysis.
        male_fhet_thresh (float): The minimum Fstat threshold for calling a sample "Male".
        female_fhet_thresh (float): The maximum Fstat threshold for calling a sample "Female".
        maf (Optional[hl.NumericExpression]): A field defining the alternate allele frequency for each row. If
        ``None``, AAF will be computed from `call`.
        reported_is_female (Optional[hl.BooleanExpression]): A field from a dataset where ``True`` means the sample is "Female" and ``False`` means "Male".

    Returns:
        A Hail Table with sex information annotated.
    """
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    mt = hl.split_multi(mt)
    mt = mt.annotate_entries(GT=hl.if_else(mt.GT.is_haploid(), hl.call(mt.GT[0], mt.GT[0]), mt.GT))

    imputed_sex = hl.impute_sex(mt.GT,
                                include_par=False,
                                aaf_threshold=maf_threshold,
                                aaf=maf,
                                female_threshold=female_fhet_thresh,
                                male_threshold=male_fhet_thresh)

    if reported_is_female is not None:
        reported_is_female_key = list(reported_is_female._indices.source.key)
        name, reported_is_female = reported_is_female._to_relational_preserving_rows_and_cols('reported_is_female')
        reported_is_female = reported_is_female.key_by(*reported_is_female_key)
        imputed_sex = imputed_sex.annotate(reported_is_female=reported_is_female[imputed_sex.key][name])
    else:
        imputed_sex = imputed_sex.annotate(reported_is_female=hl.null(hl.tbool))

    sex_check = hl.or_missing(hl.is_defined(imputed_sex.reported_is_female) & hl.is_defined(imputed_sex.is_female),
                              imputed_sex.is_female == imputed_sex.reported_is_female)

    is_passing = hl.if_else(hl.is_defined(imputed_sex.reported_is_female) & hl.is_defined(imputed_sex.is_female),
                            imputed_sex.is_female == imputed_sex.reported_is_female,
                            True)

    imputed_sex = imputed_sex.select(sex_info=hl.struct(is_female=imputed_sex.is_female,
                                                        is_passing=is_passing,
                                                        sex_check=sex_check,
                                                        f_stat=imputed_sex.f_stat,
                                                        n_called=imputed_sex.n_called,
                                                        expected_homs=imputed_sex.expected_homs,
                                                        observed_homs=imputed_sex.observed_homs))

    return imputed_sex


def identify_outliers_in_qc_stats(*,
                                  mt: hl.MatrixTable,
                                  ancestry_pop: Optional[hl.StringExpression] = None,
                                  ) -> hl.Table:
    """Identify outliers in QC statistics.

    All statistics computed with ``hl.sample_qc``.

    The result is a Hail Table with the following fields:
      * mt.col_key: The column key of the input MatrixTable
      * sample_qc_metrics:
        * is_passing: Passes every statistic.
        * r_ti_tv: Transition / Transversion ratio.
        * n_singleton: Number of singletons.
        * n_insertion: Number of insertions.
        * n_deletion: Number of deletions.
        * n_transition: Number of transitions.
        * n_transversion: Number of transversions.
        * r_het_hom_var: Ratio of heterozygotes to number of homozygote variants.
        * r_insertion_deletion: Ratio of insertions to deletions.
        * fail_r_ti_tv: An outlier in ``r_ti_tv``.
        * fail_n_singleton: An outlier in ``n_singleton``.
        * fail_n_insertion: An outlier in ``n_insertion``.
        * fail_n_deletion: An outlier in ``n_deletion``.
        * fail_n_transition: An outlier in ``n_transition``.
        * fail_n_transversion: An outlier in ``n_transversion``.
        * fail_r_het_hom_var: An outlier in ``r_het_hom_var``.
        * fail_r_insertion_deletion: An outlier in ``r_insertion_deletion``.

    Args:
        mt (MatrixTable): A Hail MatrixTable to compute sample QC statistics on. This should be LD-pruned, high quality variants.
        ancestry_pop (Optional[StringExpression]): A Hail Expression with the ancestry population label for each sample for computing outliers based on ancestry-stratified statistics.

    Returns:
        A Hail Table with sample QC statistics as well as outliers flagged for each metric.
    """
    t = hl.sample_qc(mt).cols().flatten()
    t = t.key_by(*list(mt.col_key))

    if ancestry_pop is not None:
        ancestry_pop_key = list(ancestry_pop._indices.source.key)
        name, ancestry = ancestry_pop._to_relational_preserving_rows_and_cols('pop_name')
        ancestry = ancestry.key_by(*ancestry_pop_key)
        t = t.annotate(qc_pop=ancestry[t.key][name])
    else:
        t = t.annotate(qc_pop=hl.null(hl.tstr))

    stratified_metrics = compute_stratified_metrics_filter(
        t,
        qc_metrics={metric: t[f'sample_qc.{metric}'] for metric in filtering_qc_metrics},
        strata={"qc_pop": t.qc_pop},
        metric_threshold={'sample_qc.n_singleton': (4.0, 8.0)}
        ## we'll change the singleton upper bound to 8, to better fit a 0-bounded distribution
    )

    is_passing = hl.all(*[~stratified_metrics[f'fail_{metric}'] for metric in filtering_qc_metrics])

    stratified_metrics = stratified_metrics.select(sample_qc_metrics=hl.struct(
        is_passing=is_passing,
        fail_r_ti_tv=stratified_metrics['fail_r_ti_tv'],
        fail_n_singleton=stratified_metrics['fail_n_singleton'],
        fail_n_insertion=stratified_metrics['fail_n_insertion'],
        fail_n_deletion=stratified_metrics['fail_n_deletion'],
        fail_n_transition=stratified_metrics['fail_n_transition'],
        fail_n_transversion=stratified_metrics['fail_n_transversion'],
        fail_r_het_hom_var=stratified_metrics['fail_r_het_hom_var'],
        fail_r_insertion_deletion=stratified_metrics['fail_r_insertion_deletion']
    ))

    stratified_metrics = stratified_metrics.join(t, how='outer')

    stratified_metrics = stratified_metrics.select(sample_qc_metrics=stratified_metrics.sample_qc_metrics.annotate(
        r_ti_tv=stratified_metrics['sample_qc.r_ti_tv'],
        n_singleton=stratified_metrics['sample_qc.n_singleton'],
        n_insertion=stratified_metrics['sample_qc.n_insertion'],
        n_deletion=stratified_metrics['sample_qc.n_deletion'],
        n_transition=stratified_metrics['sample_qc.n_transition'],
        n_transversion=stratified_metrics['sample_qc.n_transversion'],
        r_het_hom_var=stratified_metrics['sample_qc.r_het_hom_var'],
        r_insertion_deletion=stratified_metrics['sample_qc.r_insertion_deletion']
    ))

    return stratified_metrics


def select_high_quality_common_sites(dataset: Union[hl.MatrixTable, hl.vds.VariantDataset],
                                     is_gatk: bool,
                                     *,
                                     low_complexity_regions: hl.Table,
                                     exome_regions: hl.Table,
                                     passing_variants: hl.Table,
                                     dp_thresh: int = 10,
                                     gq_thresh: int = 20,
                                     ab_thresh: float = 0.2,
                                     maf_thresh: float = 0.01,
                                     mac_thresh: int = 10,
                                     call_rate_thresh: float = 0.95) -> hl.Table:
    """Select high quality sites for downstream sample QC.

    Args:
        dataset (MatrixTable): The MatrixTable to select high quality sites from.
        is_gatk (bool): ``True`` if the data was called using GATK. ``False`` if the data was called using DRAGEN.
        low_complexity_regions (Table): A Hail Table with intervals representing regions of the genome to filter out.
        exome_regions (Table): A Hail Table with intervals representing regions of the genome to keep.
        passing_variants (Table): A Hail Table with VQSR passing variants (or based on filters) to keep.
        dp_thresh (int): The minimum DP threshold for which a genotype is kept. DP is the sum of AD if ``is_gatk=False``.
        gq_thresh (int): The minimum GQ threshold for a kept genotype.
        ab_thresh (float): The allowable allelic balance range for a kept heterozygote genotype.
        maf_thresh (float): The minimum minor allele frequency as computed by `hl.variant_qc` in order for a variant to be considered to be common.
        mac_thresh (int): The minimum number of minor alleles observed as computed by `hl.variant_qc` in order for a variant to be considered common.
        call_rate_thresh (float): The minimum call rate of a variant as computed by `hl.variant_qc` in order for a variant to be considered.

    Returns:
          A Table with locus and alleles for common, LD-pruned high quality sites.
    """
    dataset = filter_to_high_quality_exome_intervals(dataset,
                                                     low_complexity_regions=low_complexity_regions,
                                                     exome_regions=exome_regions,
                                                     passing_variants=passing_variants,
                                                     is_gatk=is_gatk)

    if isinstance(dataset, hl.vds.VariantDataset):
        variant_data = dataset.variant_data

        if 'GT' not in variant_data.entry and 'LGT' in variant_data.entry and 'LA' in variant_data.entry:
            variant_data = variant_data.annotate_entries(
                GT=hl.vds.lgt_to_gt(variant_data.LGT, variant_data.LA)
            )
            
        if 'AD' not in variant_data.entry and 'LAD' in variant_data.entry and 'LA' in variant_data.entry:
            variant_data = variant_data.annotate_entries(AD=hl.vds.local_to_global(
                variant_data.LAD,
                variant_data.LA,
                n_alleles=hl.len(variant_data.alleles),
                fill_value=0,
                number='R'))   

        variant_data = variant_data.annotate_entries(GT=hl.coalesce(variant_data.GT, hl.call(0, 0)),
                                                     GT_backup=variant_data.GT)

        variant_data = filter_to_high_quality_common_sites(variant_data,
                                                           maf_thresh=0,
                                                           mac_thresh=mac_thresh,
                                                           call_rate_thresh=call_rate_thresh)

        variant_data = variant_data.transmute_entries(GT=variant_data.GT_backup)

        dataset = hl.vds.VariantDataset(dataset.reference_data, variant_data)

        dataset = hl.vds.to_dense_mt(dataset)

    dataset = filter_genotypes(dataset, is_gatk=is_gatk, dp_thresh=dp_thresh, gq_thresh=gq_thresh, ab_thresh=ab_thresh)

    dataset = filter_to_high_quality_common_sites(dataset,
                                                  maf_thresh=maf_thresh,
                                                  mac_thresh=mac_thresh,
                                                  call_rate_thresh=call_rate_thresh)

    return dataset.rows().select()


def infer_ancestry(dataset: Union[hl.MatrixTable, hl.vds.VariantDataset],
                   *,
                   is_gatk: bool = True,
                   dp_thresh: int = 10,
                   gq_thresh: int = 20,
                   ab_thresh: float = 0.2,
                   pca_loadings: str = 'gs://gcp-public-data--gnomad/release/4.0/pca/gnomad.v4.0.pca_loadings.ht',
                   onnx_rf: str = 'gs://gcp-public-data--gnomad/release/4.0/pca/gnomad.v4.0.RF_fit.onnx',
                   requester_pays_config: Optional[GCSRequesterPaysConfiguration] = None,
                   num_pcs: int = 20,
                   min_prob: float = 0.75,
                   log: Optional[logging.Logger] = None
                   ) -> hl.Table:
    """Use a Random Forest Model derived from gnomAD data to infer ancestry population labels.

    Does basic QC automatically:
    - Removes genotypes that do not meet the DP or GQ threshold
    - Removes hetrozygote genotypes that do not meet the AB threshold
    - Removes variants with low call rates

    The result is a Hail Table with the following fields:
      * mt.col_key: The column key of the input MatrixTable
      * ancestry:
        * ancestry_pop: The ancestry population label.

    Args:
        dataset (MatrixTable, VariantDataset): The MatrixTable with the samples to infer ancestry from. This dataset should be filtered to high quality variants.
        is_gatk (bool): True if data was called with GATK.
        dp_thresh (int): The minimum DP threshold for which a genotype is kept. DP is the sum of AD if ``is_gatk=False``.
        gq_thresh (int): The minimum GQ threshold for a kept genotype.
        ab_thresh (float): The allowable allelic balance range for a kept heterozygote genotype.
        pca_loadings (str): Path to the gnomAD PCA loadings Hail Table.
        onnx_rf (str): Path to the gnomAD random forest model.
        num_pcs (int): Number of PCs to use when applying the random forest model.
        min_prob (float): The minimum probability for a predicted population before giving a sample that label.
        log (Optional[logging.Logger]): An optional logging object.

    Returns:
        A Hail Table with two columns: sample ID and a nested struct called ``ancestry`` with one field ``ancestry_pop``.
    """
    if log is None:
        log = logging.getLogger()

    with hfs.open(onnx_rf, "rb", requester_pays_config=requester_pays_config) as f:
        data = f.read()
        onnx_fit = onnx.load_model_from_string(data)

    loadings_ht = hl.read_table(pca_loadings)

    if isinstance(dataset, hl.MatrixTable):
        dataset = dataset.semi_join_rows(loadings_ht)
    else:
        assert isinstance(dataset, hl.vds.VariantDataset)
        variant_data = dataset.variant_data
        variant_data = variant_data.semi_join_rows(loadings_ht)
        # dataset = hl.vds.filter_variants(loadings_ht)
        
        dataset = hl.vds.VariantDataset(dataset.reference_data, variant_data)
        dataset = hl.vds.to_dense_mt(dataset)

    if 'GT' not in dataset.entry and 'LGT' in dataset.entry and 'LA' in dataset.entry:
        dataset = dataset.annotate_entries(
            GT=hl.vds.lgt_to_gt(dataset.LGT, dataset.LA)
        )

    if 'AD' not in dataset.entry and 'LAD' in dataset.entry and 'LA' in dataset.entry:
        dataset = dataset.annotate_entries(AD=hl.vds.local_to_global(
            dataset.LAD,
            dataset.LA,
            n_alleles=hl.len(dataset.alleles),
            fill_value=0,
            number='R'))

    dataset = filter_genotypes(dataset, is_gatk=is_gatk, ab_thresh=ab_thresh, dp_thresh=dp_thresh, gq_thresh=gq_thresh)
    dataset = filter_to_high_quality_common_sites(dataset, maf_thresh=0, mac_thresh=0)

    dataset = dataset.select_cols().select_rows().select_entries('GT')

    dataset = dataset.persist()

    n_variants = dataset.count_rows()
    n_loadings = loadings_ht.count()

    log.info(f'Using {n_variants} variants for PC projection from {n_loadings} possible loading variants.')

    # Project new genotypes onto loadings.
    pcs_ht = hl.experimental.pc_project(
        dataset.GT, loadings_ht.loadings, loadings_ht.pca_af,
    )

    n_missing = pcs_ht.filter(hl.is_missing(pcs_ht.scores)).count()
    if n_missing > 0:
        log.info(f'missing scores for {n_missing} samples!')

    pcs_ht = pcs_ht.filter(hl.is_defined(pcs_ht.scores))

    ancestry_pops, _ = assign_population_pcs(
        pcs_ht,
        pc_cols=pcs_ht.scores[:num_pcs],
        fit=onnx_fit,
        min_prob=min_prob,
        apply_model_func=apply_onnx_classification_model,
        output_col='ancestry_pop',
    )

    ancestry_pops = ancestry_pops.select(ancestry=hl.struct(ancestry_pop=ancestry_pops.ancestry_pop))

    ancestry_pops = dataset.cols().join(ancestry_pops, how='left')

    dataset.unpersist()

    return ancestry_pops.annotate(ancestry=hl.struct(ancestry_pop=hl.coalesce(ancestry_pops.ancestry.ancestry_pop, 'unknown')))


def prune_variants(*,
                   mt: hl.MatrixTable,
                   pre_pruned_variants: Optional[hl.Table] = None,
                   r2_thresh: float = 0.2,
                   bp_window_size: int = 1_000_000) -> hl.MatrixTable:
    """LD Prune variants.

    Args:
        mt (MatrixTable): A Hail MatrixTable to prune.
        pre_pruned_variants (Optional(Table)): A list of pre-pruned variants.
        r2_thresh (float): The R2 threshold to use when pruning variants.
        bp_window_size (int): The window for considering variants to prune in base pairs.

    Returns:
        A Hail MatrixTable containing data only for LD pruned variants
    """
    if pre_pruned_variants is None:
        mt = hl.split_multi(mt)
        pruned_variants = hl.ld_prune(mt.GT, r2=r2_thresh, bp_window_size=bp_window_size)
        pruned_variants = pruned_variants.persist()
    else:
        pruned_variants = pre_pruned_variants

    mt = mt.semi_join_rows(pruned_variants)

    return mt.persist()


def select_passing_variants_from_filters(dataset: Union[hl.MatrixTable, hl.vds.VariantDataset]) -> hl.Table:
    """Select passing variants based on filters.

    Notes:

    This function supports the following cases:
    - If the filters is a boolean, treat as "PASS"
    - If the filters is "PASS", treat as "PASS"
    - If the filters is an empty array or set, treat as "PASS"

    Args:
        dataset (MatrixTable, VariantDataset) A Hail Matrix Table or VariantDataset to select passing variants.

    Returns:
        A Hail Table with just passing variants and no row annotations.
    """

    if isinstance(dataset, hl.vds.VariantDataset):
        dataset = dataset.variant_data

    if 'filters' not in dataset.row:
        dataset = dataset.annotate_rows(filters=hl.empty_set(hl.tstr))

    if dataset.row.filters.dtype != hl.tset(hl.tstr):
        dataset = dataset.annotate_rows(
            filters=hl.bind(
                lambda f: hl.case()
                .when(hl.is_missing(f), hl.empty_set(hl.tstr))
                .when(hl.is_defined(f) & (f.dtype == hl.tbool),
                      hl.empty_set(hl.tstr))  # treat any bool as PASS
                .when(hl.is_defined(f) & (f.dtype == hl.tstr),
                      hl.cond(f == 'PASS', hl.empty_set(hl.tstr), hl.set([f])))
                .when(hl.is_defined(f) & (f.dtype == hl.tarray(hl.tstr)),
                      hl.set(f))
                .default(hl.empty_set(hl.tstr)),
                dataset.filters
            )
        )

    dataset = dataset.filter_rows(dataset.filters == hl.empty_set(hl.tstr), keep=True)

    return dataset.select_rows().rows()


def calculate_sample_qc_stats(*,
                              dataset: Union[hl.MatrixTable, hl.vds.VariantDataset],
                              high_quality_sites: hl.Table,
                              pre_pruned_variants: Optional[hl.Table] = None,
                              ld_prune_r2_thresh: Optional[float] = 0.2,
                              ld_prune_bp_window_size: Optional[int] = 1_000_000,
                              chimera_rate: Optional[hl.NumericExpression] = None,
                              chimera_threshold: float = 0.05,
                              contamination_charr_thresh: float = 0.05,
                              contamination_min_af: float = 0.05,
                              contamination_max_af: float = 0.95,
                              contamination_min_dp: int = 10,
                              contamination_max_dp: int = 100,
                              contamination_min_gq: int = 20,
                              contamination_ref_AF: Optional[hl.Float64Expression] = None,
                              pcs_k: int = 10,
                              ancestry_pop_labels: Optional[hl.Table] = None,
                              # relatedness_min_individual_maf: float = 0.01,
                              # relatedness_k: Optional[int] = 10,
                              # relatedness_scores_expr: Optional[hl.ArrayNumericExpression] = None,
                              # relatedness_min_kinship: float = 0.125,
                              coverage_gq_thresh: int = 20,
                              coverage_fraction_thresh: float = 0.9,
                              sex_reported_is_female: Optional[hl.BooleanExpression] = None,
                              sex_male_fhet_thresh: float = 0.8,
                              sex_female_fhet_thresh: float = 0.2,
                              log: Optional[logging.Logger] = None) -> 'SampleQCResult':
    """Calculate Sample QC statistics from a Hail MatrixTable.

    Args:
        dataset (MatrixTable, VariantDataset): A Hail MatrixTable or VariantDataset to compute sample QC statistics on.
        high_quality_sites (Table): A Hail Table with the list of high quality sites (common variants, in exome regions, not in low complexity regions).
        pre_pruned_variants (Optional[Table]): The path to a table with LD pruned variants.
        ld_prune_r2_thresh (float): The LD-Prune R2 value when filtering out variants in linkage disequilibrium.
        ld_prune_bp_window_size (int): The window size for LD pruning variants.
        chimera_rate (NumericExpression): A numeric expression containing the chimeric read rate.
        chimera_threshold (float): The maximum threshold for which to consider a sample passing.
        contamination_charr_thresh (float): The threshold for which a sample is considered to be passing.
        contamination_min_af (float): The minimum allele frequency when using the ``compute_charr`` function in Hail.
        contamination_max_af (float): The maximum allele frequency when using the ``compute_charr`` function in Hail.
        contamination_min_dp (int): The minimum depth when using the ``compute_charr`` function in Hail.
        contamination_max_dp (int): The maximum depth when using the ``compute_charr`` function in Hail.
        contamination_min_gq (int): The minimum GQ when using the ``compute_charr`` function in Hail.
        contamination_ref_AF (Float64Expression): A float row field on the MatrixTable with the reference allele frequency when using the ``compute_charr`` function in Hail.
        pcs_k (int): The number of principal components to compute.
        ancestry_pop_labels (Optional(Table)): A Hail Table with ancestry population labels such as projected from gnomAD. There should be a column ``ancestry`` with one field ``ancestry_pop``.
        coverage_gq_thresh (int): The minimum GQ for a genotype to be considered passing.
        coverage_fraction_thresh (float): The minimum rate at which GQ is above the minimum threshold.
        sex_reported_is_female (Optional[hl.BooleanExpression]): A field from a dataset where ``True`` means the sample is "Female" and ``False`` means "Male".
        sex_male_fhet_thresh (float): The minimum Fstat threshold for calling a sample "Male".
        sex_female_fhet_thresh (float): The maximum Fstat threshold for calling a sample "Female".
        log (Optional[logging.Logger]): An optional logging object.

    Returns:
        A :class:`.SampleQCResult` with plotting and exporting functionality.
    """

    if log is None:
        log = logging.getLogger()

    log.info('subsetting the data to high quality sites')

    if isinstance(dataset, hl.MatrixTable):
        if 'DP' not in dataset.entry:
            dataset = dataset.annotate_entries(DP=hl.sum(dataset.AD))

        hq_mt = dataset.semi_join_rows(high_quality_sites).select_rows().select_entries('GT', 'AD', 'GQ', 'DP')
        mt = dataset
    else:
        assert isinstance(dataset, hl.vds.VariantDataset)
        variant_data = dataset.variant_data
        variant_data = variant_data.semi_join_rows(high_quality_sites)
        hq_vds = hl.vds.VariantDataset(dataset.reference_data, variant_data)
        # hq_vds = hl.vds.filter_variants(dataset, high_quality_sites, keep=True)
        hq_vds_variant_data = hq_vds.variant_data

        hq_vds = hl.vds.VariantDataset(hq_vds.reference_data, hq_vds_variant_data)

        hq_mt = hl.vds.to_dense_mt(hq_vds)

        if 'GT' not in hq_mt.entry and 'LGT' in hq_mt.entry and 'LA' in hq_mt.entry:
            hq_mt = hq_mt.annotate_entries(
                GT=hl.vds.lgt_to_gt(hq_mt.LGT, hq_mt.LA)
            )

        if 'AD' not in hq_mt.entry and 'LAD' in hq_mt.entry and 'LA' in hq_mt.entry:
            hq_mt = hq_mt.annotate_entries(AD=hl.vds.local_to_global(
                hq_mt.LAD,
                hq_mt.LA,
                n_alleles=hl.len(hq_mt.alleles),
                fill_value=0,
                number='R'))

        if 'DP' not in hq_mt.entry:
            hq_mt = hq_mt.annotate_entries(DP=hl.sum(hq_mt.AD))
            
        hq_mt = hq_mt.select_entries('GT', 'AD', 'GQ', 'DP')
        hq_mt = hq_mt.select_rows()

        mt = variant_data

    hq_mt = hq_mt.persist()

    autosomal_hq_mt = hq_mt.filter_rows(hq_mt.locus.in_autosome())
    chrx_hq_mt = hq_mt.filter_rows(hq_mt.locus.in_x_nonpar())

    log.info('LD pruning variants')

    pruned_hq_mt = prune_variants(mt=autosomal_hq_mt,
                                  pre_pruned_variants=pre_pruned_variants,
                                  r2_thresh=ld_prune_r2_thresh,
                                  bp_window_size=ld_prune_bp_window_size)

    n_variants, n_samples = hq_mt.count()
    n_chrx_variants = chrx_hq_mt.count_rows()
    n_pruned_variants = pruned_hq_mt.count_rows()

    samples = hq_mt.cols().select().persist()

    log.info(f'Running QC from {n_variants} autosomal variants, {n_chrx_variants} chrX variants, {n_pruned_variants} pruned variants, and {n_samples} samples.')

    if chimera_rate is not None:
        log.info('filtering chimeric reads.')
        chimera_rates = filter_chimeric_reads(chimera_rate=chimera_rate, threshold=chimera_threshold)
        chimera_rates = chimera_rates.semi_join(samples)
    else:
        chimera_rates = hq_mt.cols().select(chimera_reads=hl.struct(chimera_rate=hl.null(hl.tfloat64),
                                                                    is_passing=hl.null(hl.tbool)))

    if contamination_ref_AF is None and n_samples < 10_000:
        gnomad_sites = hl.experimental.load_dataset('gnomad_genome_sites',
                                                    version='3.1.2',
                                                    reference_genome='GRCh38',
                                                    region='us-central1',
                                                    cloud = 'gcp')

        contamination_ref_AF = 1 - gnomad_sites[autosomal_hq_mt.row_key].freq[0].AF

    log.info('using CHARR to estimate contamination rates.')

    contamination_rates = filter_contamination_rate(mt=autosomal_hq_mt,
                                                    charr_thresh=contamination_charr_thresh,
                                                    min_af=contamination_min_af,
                                                    max_af=contamination_max_af,
                                                    min_dp=contamination_min_dp,
                                                    max_dp=contamination_max_dp,
                                                    min_gq=contamination_min_gq,
                                                    ref_AF=contamination_ref_AF)

    log.info('computing principal components.')

    pcs = compute_pcs(pruned_hq_mt, k=pcs_k)

    # log.info('estimating relatedness.')
    #
    # relatedness = infer_relatedness(pruned_hq_mt,
    #                                 min_individual_maf=relatedness_min_individual_maf,
    #                                 k=relatedness_k,
    #                                 scores_expr=relatedness_scores_expr,
    #                                 min_kinship=relatedness_min_kinship)

    log.info('estimating call rate via GQ statistics.')

    call_rate_via_gq_stats = estimate_call_rate_through_gq_stats(autosomal_hq_mt,
                                                                 gq_thresh=coverage_gq_thresh,
                                                                 fraction_thresh=coverage_fraction_thresh)

    log.info('identifying outliers in QC statistics')

    if ancestry_pop_labels is not None:
        outliers = identify_outliers_in_qc_stats(mt=mt, ancestry_pop=ancestry_pop_labels.ancestry.ancestry_pop)
    else:
        outliers = identify_outliers_in_qc_stats(mt=mt)

    log.info('imputing and checking sex.')

    sex_info = impute_and_check_sex(chrx_hq_mt,
                                    male_fhet_thresh=sex_male_fhet_thresh,
                                    female_fhet_thresh=sex_female_fhet_thresh,
                                    reported_is_female=sex_reported_is_female)
    sex_info = sex_info.semi_join(samples)

    result_tables = [chimera_rates,
                     pcs,
                     contamination_rates,
                     ancestry_pop_labels,
                     # relatedness,
                     sex_info,
                     call_rate_via_gq_stats,
                     outliers
                     ]

    result = result_tables[0]
    for result_t in result_tables[1:]:
        result = result.join(result_t, how='outer')

    hq_mt.unpersist()
    pruned_hq_mt.unpersist()
    samples.unpersist()

    return SampleQCResult._from_table(result)


class SampleQCResult:
    """An interface for working with Sample QC Results."""

    @staticmethod
    def load(path: str) -> 'SampleQCResult':
        """Load a :class:`.SampleQCResult` from a Hail Table on disk.

        Args:
            path (str): The path at which the Hail Table resides.

        Returns:
            A :class:`.SampleQCResult` object.
        """
        t = hl.read_table(path)
        return SampleQCResult(t)

    @staticmethod
    def _from_table(t: hl.Table):
        t = t.persist()
        return SampleQCResult(t)

    def __init__(self, table: hl.Table):
        self._table = table

    @property
    def table(self):
        """Underlying Hail Table containing the Sample QC statistics.

        Examples:

        >>> t = qc_result.table
        >>> t = t.annotate(my_field=...)

        Returns:
            A Hail Table with Sample QC statistics.
        """
        return self._table

    def close(self):
        """Unpersist any temporary tables."""
        self.table.unpersist()

    def write(self, path: str, *, overwrite: bool = True):
        """Write the underlying table to disk.

        Args:
            path (str): The path to write the Hail Table to.
            overwrite (bool): Whether to overwrite the file if it already exists.
        """
        self.table.write(path, overwrite=overwrite)

    @property
    def contamination(self) -> hl.StructExpression:
        """A HailExpression that contains the contamination results.

        The result is a StructExpression with the following fields:

            * charr: The CHARR statistic.
            * is_passing: ``charr`` is less than the ``charr_thresh``

        See :py:func:`filter_contamination_rate` for the implementation and more details.

        Returns:
            A Hail StructExpression with the schema defined above.
        """
        return self.table.contamination

    # @property
    # def relatedness(self) -> hl.StructExpression:
    #     """A HailExpression that contains the relatedness results.
    #
    #     The result is a StructExpression with the following fields:
    #
    #         * is_passing: An unrelated sample for downstream analysis.
    #
    #
    #     See :py:func:`infer_relatedness` for the implementation and more details.
    #
    #     Returns:
    #         A Hail StructExpression with the schema defined above.
    #     """
    #     return self.table.relatedness

    @property
    def chimera_reads(self) -> hl.StructExpression:
        """A HailExpression that contains the chimeric rate results.

        The result is a StructExpression with the following fields:

            * chimera_rate: The rate of chimera reads.
            * is_passing: The rate of chimera reads is below ``threshold``.

        See :py:func:`filter_chimeric_reads` for the implementation and more details.

        Returns:
            A Hail StructExpression with the schema defined above.
        """
        return self.table.chimera_reads

    @property
    def qc_metrics(self) -> hl.StructExpression:
        """A HailExpression that contains the sample qc metrics results with outliers flagged.

        The result is a StructExpression with the following fields:

            * is_passing: Passes every statistic.
            * r_ti_tv: Transition / Transversion ratio.
            * n_singleton: Number of singletons.
            * n_insertion: Number of insertions.
            * n_deletion: Number of deletions.
            * n_transition: Number of transitions.
            * n_transversion: Number of transversions.
            * r_het_hom_var: Ratio of heterozygotes to number of homozygote variants.
            * r_insertion_deletion: Ratio of insertions to deletions.
            * fail_r_ti_tv: An outlier in ``r_ti_tv``.
            * fail_n_singleton: An outlier in ``n_singleton``.
            * fail_n_insertion: An outlier in ``n_insertion``.
            * fail_n_deletion: An outlier in ``n_deletion``.
            * fail_n_transition: An outlier in ``n_transition``.
            * fail_n_transversion: An outlier in ``n_transversion``.
            * fail_r_het_hom_var: An outlier in ``r_het_hom_var``.
            * fail_r_insertion_deletion: An outlier in ``r_insertion_deletion``.

        See :py:func:`identify_outliers_in_qc_stats` for the implementation and more details.

        Returns:
            A Hail StructExpression with the schema defined above.
        """
        return self.table.sample_qc_metrics

    @property
    def gq_fraction(self) -> hl.StructExpression:
        """A HailExpression that contains the GQ fraction results as a proxy for call rate.

        The result is a StructExpression with the following fields:

            * gq_fraction: The percentage of genotypes with GQ >= ``gq_thresh``.
            * is_passing: ``gq_fraction`` is greater than the ``fraction_thresh``

        See :py:func:`estimate_call_rate_through_gq_stats` for the implementation and more details.

        Returns:
            A Hail StructExpression with the schema defined above.
        """
        return self.table.gq_fraction

    @property
    def sex_info(self) -> hl.StructExpression:
        """A HailExpression that contains the imputed sex and sex check results.

        The result is a StructExpression with the following fields:

            * is_female: The imputed sex. ``True`` is for females, ``False`` is for males.
            * is_passing: Either ``reported sex == imputed sex`` or ``True``.
            * sex_check: ``reported sex == imputed sex`` or ``Null``.
            * f_stat: Inbreeding coefficient on the non-PAR X chromosome.
            * n_called: Number of genotypes considered.
            * expected_homs: Expected number of homozygotes.
            * observed_homs: Observed number of heterozygotes.

        See :py:func:`impute_and_check_sex` for the implementation and more details.

        Returns:
            A Hail StructExpression with the schema defined above.
        """
        return self.table.sex_info

    @property
    def ancestry_pop(self) -> hl.StringExpression:
        """A HailExpression that contains the imputed ancestry population.

        The result is a StringExpression.

        See :py:func:`infer_ancestry` for the implementation and more details.

        Returns:
            A Hail StringExpression with the ancestry population label.
        """
        return self.table.ancestry.ancestry_pop

    @property
    def pcs(self) -> hl.StructExpression:
        """A HailExpression that contains the principal components.

        The result is a StructExpression with the following fields:

            * scores: An array with the top ``k`` principal components.

        See :py:func:`compute_pcs` for the implementation and more details.

        Returns:
            A Hail StructExpression with the schema defined above.
        """
        return self.table.pcs

    @property
    def sample_id(self) -> hl.StructExpression:
        """A HailExpression that contains the key of the sample qc results.

        The key in most cases is a single field `s` with the sample ID.

        Returns:
            A Hail StructExpression with the key of the sample QC results.
        """
        return self.table.key

    def export(self,
               path: str,
               exprs: List[Union[str, hl.Expression]],
               named_exprs: Dict[str, hl.Expression],
               *,
               delimiter: str = "\t",
               missing: str = "NA",
               header: bool = True):
        """Export the sample QC results to a file.

        Examples:

            >>> sample_qc_results.export('gs://my-bucket/analysis/sample_qc.tsv',
            ...                          exprs=[sample_qc_results.relatedness.is_passing],
            ...                          named_exprs={'r_ti_tv': sample_qc_results.qc_metrics.r_ti_tv,
            ...                                       'PC1': sample_qc_results.pcs.scores[0]})

        Args:
            path (str): The file to export data to.
            exprs (List[hl.Expression]): A list of expressions to export.
            named_exprs (Dict[str, hl.Expression]): Named expressions to export.
            delimiter (str): The delimiter of the resulting exported file.
            missing (str): The value to use for missing values.
            header (bool): Include a header line.
        """
        t = self.table.select(*exprs, **named_exprs).flatten()
        t.export(path, delimiter=delimiter, missing=missing, header=header)

    def passing_samples(self, is_passing: Optional[hl.BooleanExpression] = None) -> 'SampleQCResult':
        """Filter the Sample QC results to passing samples.

        The default implementation for defining a passing sample is if it meets the following criteria:
            1. ``contamination.is_passing`` is True
            2. ``chimera_reads.is_passing`` is True or is missing.
            3. ``relatedness.is_passing`` is True for unrelated samples.
            4. ``sex_info.is_passing`` is True or is missing.
            5. ``gq_fraction.is_passing`` is True.
            6. ``qc_metrics.is_passing`` is True meaning the sample is not an outlier based on ancestry stratified QC metrics.

        Args:
            is_passing (hl.BooleanExpression): A BooleanExpression using expressions derived from this :class:`.SampleQCResult`.

        Returns:
            A :class:`.SampleQCResult` with only passing samples.
        """
        if is_passing is None:
            sample_qc_is_passing = hl.all(*[~self.table.sample_qc_metrics[f'fail_{metric}']
                                            for metric in filtering_qc_metrics
                                            if 'singleton' not in metric])

            is_passing = hl.all([self.contamination.is_passing,
                                 hl.coalesce(self.chimera_reads.is_passing, True),
                                 # self.relatedness.is_passing,
                                 hl.coalesce(self.sex_info.is_passing, True),
                                 self.gq_fraction.is_passing,
                                 sample_qc_is_passing,
                                 ])

        return self.filter_to(is_passing)

    def filter_to(self, filter: Optional[hl.BooleanExpression] = None) -> 'SampleQCResult':
        """Filter the Sample QC results to a user-specified set of samples based on a BooleanExpression.

        Args:
            filter (hl.BooleanExpression): A BooleanExpression using expressions derived from this :class:`.SampleQCResult` to identify selected samples.

        Returns:
            A :class:`.SampleQCResult` with only filtered samples.
        """
        t = self.table.filter(filter)
        return SampleQCResult(t)

    def plot_pcs(self, pc_i: int, pc_j: int) -> ggplot:
        """Plot principal components.

        Args:
            pc_i (int): This is a 1-indexed principal component index. 1=First principal component.
            pc_j (int): This is a 1-indexed principal component index. 1=First principal component.

        Returns:
            A `plotnine` ggplot object.
        """
        pc_i_col = f'PC{pc_i}'
        pc_j_col = f'PC{pc_j}'

        df = self.table.select(**{pc_i_col: self.pcs.scores[pc_i - 1],
                                  pc_j_col: self.pcs.scores[pc_j - 1],
                                  'ancestry_pop': self.ancestry_pop}).to_pandas()

        data = df.dropna(subset=[pc_i_col, pc_j_col]).copy()

        # Explicitly coerce metric to float and purge pd.NA
        data[pc_i_col] = pd.to_numeric(data[pc_i_col], errors='coerce').astype(float)
        data[pc_j_col] = pd.to_numeric(data[pc_j_col], errors='coerce').astype(float)

        plot = (
                ggplot(data, aes(x=pc_i_col, y=pc_j_col, color="ancestry_pop"))
                + geom_point()
                + theme_minimal()
                + theme(legend_position="right",
                        plot_title=element_text(size=18),
                        axis_title_x=element_text(size=14),
                        axis_title_y=element_text(size=14),
                        axis_text_y=element_text(size=12),
                        legend_title=element_text(size=14),
                        legend_text=element_text(size=12),
                        panel_background=element_rect(fill='white'),
                        plot_background=element_rect(fill='white')
                        )
                + coord_fixed(ratio=1)
                + labs(title="PCA Plot", x=pc_i_col, y=pc_j_col, color='Ancestry')
        )

        return plot

    def plot_qc_metric_boxplots(self) -> Dict[str, ggplot]:
        """Plot QC metrics as boxplots.

        Possible metrics:
            * r_ti_tv: Transition / Transversion ratio.
            * n_singleton: Number of singletons.
            * n_insertion: Number of insertions.
            * n_deletion: Number of deletions.
            * n_transition: Number of transitions.
            * n_transversion: Number of transversions.
            * r_het_hom_var: Ratio of heterozygotes to number of homozygote variants.
            * r_insertion_deletion: Ratio of insertions to deletions.

        Returns:
            A Dictionary of metric name to `plotnine` ggplot object.
        """
        global filtering_qc_metrics

        df = self.table.select(self.ancestry_pop, **self.qc_metrics).to_pandas()

        df['ancestry_pop'] = df['ancestry_pop'].astype('category')

        plots = {}
        for metric in filtering_qc_metrics:
            data = df.dropna(subset=[metric]).copy()

            # Explicitly coerce metric to float and purge pd.NA
            data[metric] = pd.to_numeric(data[metric], errors='coerce').astype(float)

            data[f'fail_{metric}_str'] = data[f'fail_{metric}'].map({True: 'Removed', False: 'Retained'})

            plot = (
                    ggplot(data, aes(x="ancestry_pop", y=metric))
                    + geom_jitter(data=data[data[f'fail_{metric}']].copy(),
                                  mapping=aes(x="ancestry_pop", y=metric, color=f'fail_{metric}_str'), width=0.2, height=0, alpha=0.6, size=2)
                    + geom_boxplot(outlier_shape=None)  # hide default outliers
                    + scale_color_manual(values={'Removed': 'darkred', 'Retained': 'gray'})
                    + theme_minimal()
                    + coord_flip()
                    + theme(legend_position="right",
                            plot_title=element_text(size=18),
                            axis_title_x=element_text(size=14),
                            axis_title_y=element_text(size=14),
                            axis_text_y=element_text(size=12),
                            legend_title=element_text(size=14),
                            legend_text=element_text(size=12),
                            panel_background=element_rect(fill='white'),
                            plot_background=element_rect(fill='white')
                            )
                    + labs(title="", x="Ancestry", y=metric, color='QC Outcome'))

            plots[metric] = plot

        return plots

    def plot_qc_metric_density(self) -> Dict[str, ggplot]:
        """Plot QC metrics as density plots.

        Possible metrics:
            * r_ti_tv: Transition / Transversion ratio.
            * n_singleton: Number of singletons.
            * n_insertion: Number of insertions.
            * n_deletion: Number of deletions.
            * n_transition: Number of transitions.
            * n_transversion: Number of transversions.
            * r_het_hom_var: Ratio of heterozygotes to number of homozygote variants.
            * r_insertion_deletion: Ratio of insertions to deletions.

        Returns:
            A Dictionary of metric name to `plotnine` ggplot object.
        """
        df = self.table.select(self.ancestry_pop, **self.qc_metrics).to_pandas()

        df['ancestry_pop'] = df['ancestry_pop'].astype('category')

        plots = {}
        for metric in filtering_qc_metrics:
            data = df.dropna(subset=[metric]).copy()

            # Explicitly coerce metric to float and purge pd.NA
            data[metric] = pd.to_numeric(data[metric], errors='coerce').astype(float)

            plot = (
                    ggplot(data, aes(x=metric, color="ancestry_pop"))
                    + geom_density()
                    + theme_minimal()
                    + theme(legend_position="right",
                            plot_title=element_text(size=18),
                            axis_title_x=element_text(size=14),
                            axis_title_y=element_text(size=14),
                            axis_text_y=element_text(size=12),
                            legend_title=element_text(size=14),
                            legend_text=element_text(size=12),
                            panel_background=element_rect(fill='white'),
                            plot_background=element_rect(fill='white')
                            )
                    + labs(title="", x=metric, color='Ancestry', y='Density'))
            plots[metric] = plot

        return plots

    def plot_gq_fraction_density(self) -> ggplot:
        """Plot GQ fraction as a density plot.

        Returns:
            A `plotnine` ggplot object.
        """
        df = self.table.select(self.ancestry_pop, **self.gq_fraction).to_pandas()

        df['ancestry_pop'] = df['ancestry_pop'].astype('category')

        metric = 'gq_fraction'

        data = df.dropna(subset=[metric]).copy()

        # Explicitly coerce metric to float and purge pd.NA
        data[metric] = pd.to_numeric(data[metric], errors='coerce').astype(float)

        plot = (
                ggplot(data, aes(x=metric, color="ancestry_pop"))
                + geom_density()
                + theme_minimal()
                + theme(legend_position="right",
                        plot_title=element_text(size=18),
                        axis_title_x=element_text(size=14),
                        axis_title_y=element_text(size=14),
                        axis_text_y=element_text(size=12),
                        legend_title=element_text(size=14),
                        legend_text=element_text(size=12),
                        panel_background=element_rect(fill='white'),
                        plot_background=element_rect(fill='white')
                        )
                + labs(title="", x=metric, color='Ancestry', y='Density'))

        return plot

    def plot_chimera_reads(self) -> ggplot:
        """Plot chimera reads as a boxplot.

        Colors samples that fail the threshold as a different color.

        Returns:
            A `plotnine` ggplot object.
        """
        data = self.table.select(self.ancestry_pop, **self.chimera_reads).to_pandas()

        data['ancestry_pop'] = data['ancestry_pop'].astype('category')

        metric = 'chimera_rate'

        data = data.dropna(subset=[metric]).copy()

        # Explicitly coerce metric to float and purge pd.NA
        data[metric] = pd.to_numeric(data[metric], errors='coerce').astype(float)

        data[f'fail_{metric}_str'] = data['is_passing'].map({False: 'Removed', True: 'Retained'})

        plot = (
                ggplot(data, aes(x="ancestry_pop", y=metric))
                + geom_jitter(data[~data['is_passing']].copy(),
                              aes(width=0.2, height=0, color=f'fail_{metric}_str'), alpha=0.6,
                              size=2)
                + geom_boxplot(outlier_shape=None)  # hide default outliers
                + scale_color_manual(values={'Removed': 'darkred', 'Retained': 'gray'})
                + theme_minimal()
                + coord_flip()
                + theme(legend_position="right",
                        plot_title=element_text(size=18),
                        axis_title_x=element_text(size=14),
                        axis_title_y=element_text(size=14),
                        axis_text_y=element_text(size=12),
                        legend_title=element_text(size=14),
                        legend_text=element_text(size=12),
                        panel_background=element_rect(fill='white'),
                        plot_background=element_rect(fill='white')
                        )
                + labs(title="", x="Ancestry", y=metric, color='QC Outcome'))

        return plot
