import hail as hl

from typing import Dict, Optional


def filter_to_high_quality_exome_intervals(mt: hl.MatrixTable,
                                           low_complexity_regions: hl.Table,
                                           exome_regions: hl.Table) -> hl.MatrixTable:
    mt = hl.filter_intervals(mt, low_complexity_regions, keep=False)
    mt = hl.filter_intervals(mt, exome_regions, keep=True)
    return mt


def infer_capture_kit(mt: hl.MatrixTable,
                      capture_kit_intervals: hl.Table) -> hl.Table:
    """Do a PCA of call rates across capture kit intervals."""
    raise NotImplementedError


def filter_contamination_rate(mt: hl.MatrixTable,
                              min_af: float = 0.05,
                              max_af: float = 0.95,
                              min_dp: int = 10,
                              max_dp: int = 100,
                              min_gq: int = 20,
                              ref_AF: Optional[hl.Float64Expression] = None) -> hl.Table:
    charr_t = hl.compute_charr(mt, min_af=min_af, max_af=max_af, min_dp=min_dp, max_dp=max_dp, min_gq=min_gq, ref_AF=ref_AF)
    t = charr_t.annotate(contamination=hl.struct(charr_t.charr))
    return t


def filter_chimeric_reads(*,
                          t: hl.Table,
                          chimeric_read_rate_col: str,
                          threshold: float) -> hl.Table:
    raise NotImplementedError


def filter_genotypes(mt: hl.MatrixTable,
                     *,
                     dp_thresh: int = 10,
                     gq_thresh: int = 20,
                     ab_thresh: float = 0.2
                     ) -> hl.MatrixTable:
    mt = mt.annotate_entries(DP=hl.sum(mt.AD))

    alt_ab = mt.AD[1] / mt.DP
    is_het_ab_within_thresh = (alt_ab >= ab_thresh) & (alt_ab <= (1 - ab_thresh))

    mt = mt.filter_entries(mt.DP >= dp_thresh)
    mt = mt.filter_entries(mt.GQ >= gq_thresh)
    mt = mt.filter_entries(mt.GT.is_het() & is_het_ab_within_thresh)

    return mt


def filter_to_common_pruned_sites(mt: hl.MatrixTable,
                                  *,
                                  maf_thresh: float = 0.01,
                                  mac_thresh: int = 10,
                                  call_rate_thresh: float = 0.95,
                                  r2_thresh: float = 0.2,
                                  bp_window_size: int = 1000000) -> hl.MatrixTable:
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows((hl.min(mt.variant_qc.AF) > maf_thresh) & (mt.variant_qc.call_rate > call_rate_thresh))
    mt = mt.filter_rows(mt.variant_qc.AC >= mac_thresh)
    pruned_variants = hl.ld_prune(mt.GT, r2=r2_thresh, bp_window_size=bp_window_size)
    mt = mt.semi_join_rows(pruned_variants)
    return mt


def compute_pcs(mt: hl.MatrixTable, *, k: int = 10) -> hl.Table:
    _, scores, _ = hl.hwe_normalized_pca(mt.GT, k=k)
    scores = scores.annotate(pcs=scores.scores)
    return scores


def infer_ancestry(mt: hl.MatrixTable,
                   *,
                   reference_mts: Dict[str, hl.MatrixTable],
                   pop_labels: hl.Table,
                   call_rate_thresh: float = 0.99,
                   maf_thresh: float = 0.001,
                   r2_thresh: float = 0.1,
                   bp_window_size: int = 500000,
                   hwe_pca_normalized_k: int = 10) -> hl.Table:
    """Run PCA with reference datasets and use population labels and a random forest model to infer ancestry"""
    raise NotImplementedError


def infer_relatedness(mt: hl.MatrixTable,
                      min_individual_maf: float = 0.01,
                      k: Optional[int] = None,
                      scores_expr: Optional[hl.ArrayNumericExpression] = None,
                      min_kinship: Optional[float] = 0.125) -> hl.Table:
    pairs = hl.pc_relate(mt.GT,
                         min_individual_maf=min_individual_maf,
                         k=k,
                         scores_expr=scores_expr,
                         statistics='kin',
                         min_kinship=min_kinship)

    unrelated = hl.maximal_independent_set(pairs.i, pairs.j, keep=True)
    unrelated_samples = hl.set(unrelated.node.collect(_localize=False))

    samples = mt.cols().select()
    samples = samples.annotate(relatedness=hl.struct(include=unrelated_samples.contains(samples.key)))
    return samples


def filter_call_rate_through_gq_stats(mt: hl.MatrixTable,
                                      call_rate_thresh: float = 0.9) -> hl.Table:
    """The proportion of genotypes with calls with GQ >= 20"""
    raise NotImplementedError


def impute_and_check_sex(mt: hl.MatrixTable,
                         maf_threshold: float = 0.0,
                         male_fhet_thresh: float = 0.8,
                         female_fhet_thresh: float = 0.2,
                         maf: Optional[hl.NumericExpression] = None,
                         reported_is_female: Optional[hl.Table] = None) -> hl.Table:
    # remove PAR
    reference_genome = mt.locus.dtype.reference_genome
    par_intervals = reference_genome.par
    x_contigs = hl.set(reference_genome.x_contigs)

    mt = mt.filter_rows(x_contigs.contains(mt.locus.contig))
    mt = hl.filter_intervals(mt, par_intervals, keep=False)

    imputed_sex = hl.impute_sex(mt.GT,
                                aaf_threshold=maf_threshold,
                                aaf=maf,
                                female_threshold=female_fhet_thresh,
                                male_threshold=male_fhet_thresh)

    if reported_is_female is not None:
        imputed_sex = imputed_sex.annotate(reported_is_female=reported_is_female[imputed_sex.key].is_female)
    else:
        imputed_sex = imputed_sex.annotate(reported_is_female=hl.null(hl.tbool))

    imputed_sex = imputed_sex.select(sex_info=hl.struct(is_female=imputed_sex.is_female,
                                                        sex_check=imputed_sex.is_female == imputed_sex.reported_is_female,
                                                        f_stat=imputed_sex.f_stat,
                                                        n_called=imputed_sex.n_called,
                                                        expected_homs=imputed_sex.expected_homs,
                                                        observed_homs=imputed_sex.observed_homs))

    return imputed_sex



def compute_quality_statistics(mt: hl.MatrixTable) -> hl.MatrixTable:
    raise NotImplementedError


def select_high_quality_sites(mt: hl.MatrixTable,
                              *,
                              low_complexity_regions: hl.Table,
                              exome_regions: hl.Table,
                              dp_thresh: int = 10,
                              gq_thresh: int = 20,
                              ab_thresh: float = 0.2,
                              maf_thresh: float = 0.01,
                              mac_thresh: int = 10,
                              call_rate_thresh: float = 0.95,
                              r2_thresh: float = 0.2,
                              bp_window_size: int = 1000000) -> hl.Table:
    mt = filter_to_high_quality_exome_intervals(mt,
                                                low_complexity_regions=low_complexity_regions,
                                                exome_regions=exome_regions)

    mt = filter_genotypes(mt, dp_thresh=dp_thresh, gq_thresh=gq_thresh, ab_thresh=ab_thresh)

    mt = filter_to_common_pruned_sites(mt,
                                       maf_thresh=maf_thresh,
                                       mac_thresh=mac_thresh,
                                       call_rate_thresh=call_rate_thresh,
                                       r2_thresh=r2_thresh,
                                       bp_window_size=bp_window_size)

    return mt.rows().select()


def calculate_sample_qc_stats(mt: hl.MatrixTable,
                              high_quality_sites: hl.Table,
                              pcs_k: int = 10,
                              sex_reported_is_female: Optional[hl.Table] = None,
                              sex_male_fhet_thresh: float = 0.8,
                              sex_female_fhet_thresh: float = 0.2,
                              relatedness_min_kinship: float = 0.125,
                              ancestry_reference_mts: Optional[Dict[str, hl.MatrixTable]] = None,
                              ancestry_pop_labels: Optional[hl.Table] = None,
                              ancestry_call_rate_thresh: float = 0.99,
                              ancestry_maf_thresh: float = 0.001,
                              ancestry_r2_thresh: float = 0.1,
                              ancestry_bp_window_size: int = 500000,
                              ancestry_hwe_pca_normalized_k: int = 10
                              ) -> hl.Table:
    mt = mt.semi_join_rows(high_quality_sites)

    pcs = compute_pcs(mt, k=pcs_k)

    if ancestry_reference_mts is not None and ancestry_pop_labels is not None:
        pop_labels = infer_ancestry(mt,
                                    reference_mts=ancestry_reference_mts,
                                    pop_labels=ancestry_pop_labels,
                                    call_rate_thresh=ancestry_call_rate_thresh,
                                    maf_thresh=ancestry_maf_thresh,
                                    r2_thresh=ancestry_r2_thresh,
                                    bp_window_size=ancestry_bp_window_size,
                                    hwe_pca_normalized_k=ancestry_hwe_pca_normalized_k)
    else:
        pop_labels = mt.cols().select(ancestry=hl.null(hl.tstruct(pop=hl.tstr)))

    relatedness = infer_relatedness(mt, min_kinship=relatedness_min_kinship)

    sex_info = impute_and_check_sex(mt,
                                    male_fhet_thresh=sex_male_fhet_thresh,
                                    female_fhet_thresh=sex_female_fhet_thresh,
                                    reported_is_female=sex_reported_is_female)

    # identify outliers
    outliers = ...

    result_tables = [pcs, pop_labels, relatedness, sex_info, outliers]
    result = result_tables[0]
    for result_t in result_tables:
        result = result.join(result_t, how='outer')

    return result
