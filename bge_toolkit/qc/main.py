import logging
from typing import List, Optional, Union

import hail as hl

from ..common.io import load_data
from ..common.binned_exprs import (mac_bin, maf_bin, max_gp_bin, gq_bin, qual_approx_bin, info_score_bin, dp_bin)
from ..common.utils import apply_filters, save_ggplot, setup_logging_and_qob, setup_logger
from .joint_callset import JointCallSet
from .concordance import Statistic
from .qc import SampleQCResult, calculate_sample_qc_stats, infer_ancestry, select_high_quality_common_sites, select_passing_variants_from_filters
from .utils import ALL_AGG, ALL_GROUP, JoinType


def concordance(*,
                exome_path: str,
                imputation_path: str,
                out_dir: str,
                exome_maf: bool = False,
                exome_mac: bool = False,
                imputation_maf: bool = False,
                imputation_mac: bool = False,
                info_score: bool = False,
                dp: bool = False,
                gq: bool = False,
                max_gp: bool = False,
                qual_approx: bool = False,
                sample_list: Optional[str] = None,
                variant_list: Optional[str] = None,
                contig: Optional[List[str]] = None,
                n_samples: Optional[int] = None,
                n_variants: Optional[int] = None,
                downsample_samples: Optional[float] = None,
                downsample_variants: Optional[float] = None,
                join_type: JoinType = JoinType.inner,
                run_variants: bool = False,
                run_samples: bool = False,
                run_global: bool = False,
                log_path: Optional[str] = None,
                hail_init_kwargs: Optional[dict] = None,
                log: Optional[logging.Logger] = None,
                n_exome_partitions: Optional[int] = None,
                n_imp_partitions: Optional[int] = None):
    """Run the default implementation of concordance used in the CLI.

    Args:
        exome_path (str): Exome dataset to compare to.
        imputation_path (str): Imputation dataset to compare to.
        out_dir (str): Output directory.
        exome_maf (bool): Bin concordance counts by Exome Minor Allele Frequency.
        exome_mac (bool): Bin concordance counts by Exome Minor Allele Counts.
        imputation_maf (bool): Bin concordance counts by Imputation Minor Allele Frequency.
        imputation_mac (bool): Bin concordance counts by Imputation Minor Allele Counts.
        info_score (bool): Bin concordance counts by imputation INFO score.
        dp (bool): Bin concordance counts by exome genotype DP.
        gq (bool): Bin concordance counts by exome genotype GQ.
        max_gp (bool): Bin concordance counts by maximum value of imputation GP.
        qual_approx (bool): Bin concordance counts by variant Approx Qual Score.
        sample_list (Optional[str]): Filter to samples listed in the file.
        variant_list (Optional[str]): Filter to variants listed in the file.
        contig (List[str]): Filter to variants in the contig.
        n_samples (Optional[int]): Filter to the first N overlapping samples.
        n_variants (Optional[int]): Filter to the first N variants.
        downsample_samples (Optional[float]): Downsample to X fraction of samples.
        downsample_variants (Optional[float]): Downsample to X fraction of variants.
        join_type (JoinType): Join type.
        run_variants (bool): Run variant concordance.
        run_samples (bool): Run sample concordance.
        run_global (bool): Run global concordance.
        log_path (Optional[str]): Log file path to write to.
        hail_init_kwargs (Optional[dict]): Keyword arguments to hl.init().
        log (Optional[logging.Logger]): Logging object to log to.
        n_exome_partitions (Optional[int]): Number of partitions to load the exome dataset with if it's a MatrixTable.
        n_imp_partitions (Optional[int]): Number of partitions to load the imputation dataset with if it's a MatrixTable.
    """
    log = setup_logging_and_qob(log_name='concordance',
                                log_path=log_path,
                                hail_init_kwargs=hail_init_kwargs,
                                log=log)

    exome = load_data(path=exome_path, descriptor='exome', log=log, n_partitions=n_exome_partitions)
    imputation = load_data(path=imputation_path, descriptor='imputation', log=log, n_partitions=n_imp_partitions)

    out_dir = out_dir.rstrip('/') + '/'

    exo_row_group_by_var = {}
    imp_row_group_by_var = {}
    exo_row_agg_fields = {}
    imp_row_agg_fields = {}
    exo_entry_agg_fields = {}
    imp_entry_agg_fields = {}

    category_orders = {}

    requires_expensive_agg = False

    if exome_maf:
        if run_samples:
            log.warning(f'ignoring grouping variables and aggregation variables for sample concordance -- Exome MAF.')
        if run_variants or run_global:
            requires_expensive_agg = True
            exome = exome.annotate_rows(
                AF=hl.agg.sum(exome.GT.n_alt_alleles()) / (2 * hl.agg.count_where(hl.is_defined(exome.GT))))
            exo_row_group_by_var['exo_maf'] = maf_bin
            category_orders['exo_maf'] = ['<1%', '1-2%', '2-5%', '>5%']
            log.info('Aggregating by Exome MAF.')

    if imputation_maf:
        if run_samples:
            log.warning(f'ignoring grouping variables and aggregation variables for sample concordance -- Imputation MAF.')
        if run_variants or run_global:
            requires_expensive_agg = True
            imputation = imputation.annotate_rows(
                AF=hl.agg.sum(imputation.GT.n_alt_alleles()) / (2 * hl.agg.count_where(hl.is_defined(imputation.GT))))
            imp_row_group_by_var['imp_maf'] = maf_bin
            category_orders['imp_maf'] = ['<1%', '1-2%', '2-5%', '>5%']
            log.info('Aggregating by Imputation MAF.')

    if exome_mac:
        if run_samples:
            log.warning(f'ignoring grouping variables and aggregation variables for sample concordance -- Exome MAC.')
        if run_variants or run_global:
            requires_expensive_agg = True
            exome = exome.annotate_rows(AC=hl.agg.sum(exome.GT.n_alt_alleles()))
            exo_row_group_by_var['exo_mac'] = mac_bin
            category_orders['exo_mac'] = ['1', '2-5', '6-10', '10+']
            log.info('Aggregating by Exome MAC.')

    if imputation_mac:
        if run_samples:
            log.warning(f'ignoring grouping variables and aggregation variables for sample concordance -- Imputation MAC.')
        if run_variants or run_global:
            requires_expensive_agg = True
            imputation = imputation.annotate_rows(AC=hl.agg.sum(imputation.GT.n_alt_alleles()))
            imp_row_group_by_var['imp_mac'] = mac_bin
            category_orders['imp_mac'] = ['1', '2-5', '6-10', '10+']
            log.info('Aggregating by Imputation MAC.')

    if qual_approx:
        if run_samples:
            log.warning(f'ignoring grouping variables and aggregation variables for sample concordance -- Exome Qual Approx.')
        if run_variants or run_global:
            requires_expensive_agg = True
            assert 'QUALapprox' in exome.info
            exome = exome.annotate_rows(QUALapprox=exome.info.QUALapprox)
            exo_row_agg_fields['qual_approx'] = qual_approx_bin
            log.info('Aggregating by Exome Qual Approx.')

    if info_score:
        if run_samples:
            log.warning(f'ignoring grouping variables and aggregation variables for sample concordance -- INFO score.')
        if run_variants or run_global:
            requires_expensive_agg = True
            assert 'INFO' in imputation.info
            imputation = imputation.annotate_rows(INFO=imputation.info.INFO)
            imp_row_agg_fields['info_score'] = info_score_bin
            log.info('Aggregating by Imputation INFO score.')

    if dp:
        if run_samples:
            log.warning(f'ignoring grouping variables and aggregation variables for sample concordance -- Exome DP.')
        if run_variants or run_global:
            requires_expensive_agg = True
            if 'DP' not in exome.entry:
                exome = exome.annotate_entries(DP=hl.sum(exome.AD))
            exo_entry_agg_fields['DP'] = dp_bin
            log.info('Aggregating by Exome DP.')

    if gq:
        if run_samples:
            log.warning(f'ignoring grouping variables and aggregation variables for sample concordance -- Exome GQ.')
        if run_variants or run_global:
            requires_expensive_agg = True
            assert 'GQ' in exome.entry
            exo_entry_agg_fields['GQ'] = gq_bin
            log.info('Aggregating by Exome GQ.')

    if max_gp:
        if run_samples:
            log.warning(f'ignoring grouping variables and aggregation variables for sample concordance -- MAX GP.')
        if run_variants or run_global:
            requires_expensive_agg = True
            assert 'GP' in imputation.entry
            imputation = imputation.annotate_entries(MAX_GP=hl.max(imputation.GP))
            imp_entry_agg_fields['MAX_GP'] = max_gp_bin
            log.info('Aggregating by Imputation MAX GP.')

    exome, downsampled_exome = apply_filters(mt=exome,
                                             description='exome',
                                             log=log,
                                             sample_list=sample_list,
                                             variant_list=variant_list,
                                             contig=contig,
                                             n_samples=None,  # we want to downsample overlapping samples
                                             n_variants=n_variants,
                                             downsample_samples=None,  # we want to downsample overlapping samples
                                             downsample_variants=downsample_variants)

    imputation, downsampled_imp = apply_filters(mt=imputation,
                                                description='imputation',
                                                log=log,
                                                sample_list=sample_list,
                                                variant_list=variant_list,
                                                contig=contig,
                                                n_samples=None,  # we want to downsample overlapping samples
                                                n_variants=n_variants,
                                                downsample_samples=None,  # we want to downsample overlapping samples
                                                downsample_variants=downsample_variants)

    if requires_expensive_agg and (not downsampled_imp or not downsampled_exome):
        log.warning('Extra agregation options are expensive and no downsampling mechanism was selected.')

    log.info('Using sample join type "inner"')
    log.info(f'Using variant join type "{join_type.value}"')

    log.info('Splitting multi-allelic variants in the exome dataset.')
    exome = hl.methods.split_multi(exome)

    log.info('Splitting multi-allelic variants in the imputation dataset.')
    imputation = hl.methods.split_multi(imputation)

    if n_samples is not None or downsample_samples is not None:
        with JointCallSet(exome, imputation) as joint_callset:
            if n_samples is not None:
                log.info(f'Taking the first {n_samples} samples.')
                new_samples = joint_callset.sample_overlaps().head(n_samples)
            if downsample_samples is not None:
                log.info(f'Downsampling samples by {downsample_samples}.')
                new_samples = joint_callset.sample_overlaps().sample(n_samples, seed=10243523)

            exome = exome.semi_join_cols(new_samples)
            imputation = imputation.semi_join_cols(new_samples)

    with (JointCallSet(exome, imputation) as joint_callset):
        if run_variants:
            fields = {'concordance_rate': 'Concordance Rate',
                      'nonref_concordance_rate': 'Non-ref Concordance',
                      'n_samples': '# Samples',
                      'n_concordant': '# Concordant',
                      'n_total': '# Total Genotypes',
                      'HET_HOMREF': 'Het->Ref',
                      'HET_HET': 'Het->Het'}

            output_path = f'{out_dir}variant_overlaps.ht'
            joint_callset.variant_overlaps().write(output_path, overwrite=True)
            n_variants = joint_callset.n_variant_overlaps
            log.info(f'Found {n_variants} overlapping variants.')

            variant_results = joint_callset.variant_concordance(exo_entry_agg_fields=exo_entry_agg_fields,
                                                                imp_entry_agg_fields=imp_entry_agg_fields,
                                                                join_type=join_type)
            output_path = f'{out_dir}variant_conc.ht'
            variant_results.write(output_path=output_path, overwrite=True)
            log.info(f'Wrote the raw variant concordance data to {output_path}')

            variant_view = variant_results.get_view(agg_names=[*exo_entry_agg_fields, *imp_entry_agg_fields],
                                                    ordering=category_orders)

            variant_view.export_all(f'{out_dir}variant-results/', sep='\t', index=False, fields=fields)
            log.info(f'Exported variant results to {out_dir}variant-results/')

        if run_samples:
            fields = {'concordance_rate': 'Concordance Rate',
                      'nonref_concordance_rate': 'Non-ref Concordance',
                      'n_variants': '# Variants',
                      'n_concordant': '# Concordant',
                      'n_total': '# Total Genotypes',
                      'HET_HOMREF': 'Het->Ref',
                      'HET_HET': 'Het->Het'}

            output_path = f'{out_dir}sample_overlaps'
            sample_overlaps = joint_callset.sample_overlaps()
            sample_overlaps.write(f'{output_path}.ht', overwrite=True)
            sample_overlaps.export(f'{output_path}.tsv')
            n_samples = joint_callset.n_sample_overlaps
            log.info(f'Found {n_samples} overlapping samples.')

            sample_results = joint_callset.sample_concordance(join_type=join_type)
            output_path = f'{out_dir}sample_conc.ht'
            sample_results.write(output_path=output_path, overwrite=True)
            log.info(f'Wrote the raw sample concordance data to {output_path}')

            sample_view = sample_results.get_view()
            sample_view.export_all(f'{out_dir}sample-results/', sep='\t', index=False, fields=fields)
            log.info(f'Exported sample results to {out_dir}sample-results/')

        if run_global:
            fields = {'concordance_rate': 'Concordance Rate',
                      'nonref_concordance_rate': 'Non-ref Concordance',
                      'n_variants': '# Variants',
                      'n_samples': '# Samples',
                      'n_concordant': '# Concordant',
                      'n_total': '# Total Genotypes',
                      'HET_HOMREF': 'Het->Ref',
                      'HET_HET': 'Het->Het'}

            global_results = joint_callset.global_concordance(
                exo_row_group_by_var=exo_row_group_by_var,
                imp_row_group_by_var=imp_row_group_by_var,
                exo_row_agg_fields=exo_row_agg_fields,
                imp_row_agg_fields=imp_row_agg_fields,
                exo_entry_agg_fields=exo_entry_agg_fields,
                imp_entry_agg_fields=imp_entry_agg_fields,
                join_type=join_type
            )

            global_ht_path = f'{out_dir}global_concordance_table.ht'
            global_results.write(global_ht_path, overwrite=True)
            log.info(f'Wrote underlying concordance data to {global_ht_path}')

            global_view = global_results.get_view(group_names=[*exo_row_group_by_var, *imp_row_group_by_var],
                                                  agg_names=[*exo_row_agg_fields, *imp_row_agg_fields,
                                                             *exo_entry_agg_fields, *imp_entry_agg_fields],
                                                  ordering=category_orders)

            global_view.export_all(f'{out_dir}global-results/', sep='\t', index=False, fields=fields)
            log.info(f'Exported global results to {out_dir}global-results/')

            global_view.plot_all(out_dir=f'{out_dir}global-plots/nonref-conc/', statistic=Statistic.NONREF_CONCORDANCE)
            global_view.plot_all(out_dir=f'{out_dir}global-plots/f1-score/', statistic=Statistic.F1_SCORE)
            log.info(f'Exported global plots to {out_dir}global-plots/')

            all_view = global_results.get_view(group_names=[ALL_GROUP], agg_names=[ALL_AGG], ordering=category_orders)
            all_df = all_view.result()

            n_variants = all_df['n_variants'].tolist()[0]
            n_samples = all_df['n_samples'].tolist()[0]
            concordance_rate = all_df['concordance_rate'].tolist()[0]
            nonref_concordance_rate = all_df['nonref_concordance_rate'].tolist()[0]

            log.info(f'Found {n_variants} variants in total')
            log.info(f'Found {n_samples} samples in total')
            log.info(f'Found concordance is {float(concordance_rate):.2f}%')
            log.info(f'Found nonref concordance is {float(nonref_concordance_rate):.2f}%')


def sample_qc(*,
              dataset: Union[hl.MatrixTable, hl.vds.VariantDataset],
              out_dir: str,
              is_gatk: bool,
              exome_regions: hl.Table,
              low_complexity_regions: hl.Table,
              passing_variants: Optional[hl.Table] = None,
              hq_sites_dp_thresh: int = 10,
              hq_sites_gq_thresh: int = 20,
              hq_sites_ab_thresh: float = 0.2,
              hq_sites_maf_thresh: float = 0.01,
              hq_sites_mac_thresh: int = 10,
              hq_sites_call_rate_thresh: float = 0.95,
              pre_pruned_variants: Optional[hl.Table] = None,
              r2_thresh: float = 0.2,
              bp_window_size: int = 1000000,
              sex_reported_is_female: Optional[hl.BooleanExpression] = None,
              sex_check_male_fhet_thresh: float = 0.8,
              sex_check_female_fhet_thresh: float = 0.2,
              # relatedness_min_kinship: float = 0.125,
              pcs_k: int = 10,
              ancestry_pca_loadings: str = 'gs://gcp-public-data--gnomad/release/4.0/pca/gnomad.v4.0.pca_loadings.ht',
              ancestry_onnx_rf: str = 'gs://gcp-public-data--gnomad/release/4.0/pca/gnomad.v4.0.RF_fit.onnx',
              ancestry_num_pcs: int = 20,
              ancestry_min_prob: float = 0.75,
              chimera_rate: Optional[hl.NumericExpression] = None,
              chimera_threshold: float = 0.05,
              contamination_charr_thresh: float = 0.05,
              contamination_min_af: float = 0.05,
              contamination_max_af: float = 0.95,
              contamination_min_dp: int = 10,
              contamination_max_dp: int = 100,
              contamination_min_gq: int = 20,
              contamination_ref_af: Optional[hl.NumericExpression] = None,
              coverage_gq_thresh: int = 20,
              coverage_fraction_thresh: float = 0.9,
              log: Optional[logging.Logger] = None,
              make_plots: bool = True) -> SampleQCResult:
    """Run the default implementation of Sample QC.

    Args:
        dataset (Union(MatrixTable, VariantDataset)): Hail MatrixTable or VariantDataset to QC.
        out_dir (str): Path to output directory to write sample QC results to.
        is_gatk (bool): True if data was called with GATK.
        exome_regions (Table): Table of locus intervals specifying valid exome regions.
        low_complexity_regions (Table): Table of locus intervals specifying low complexity regions to remove.
        passing_variants (Optional[Table]): Table of passing variants from filters or VQSR.
        hq_sites_dp_thresh (int): For finding high quality sites, the minimum DP threshold for which a genotype is kept. DP is the sum of AD if ``is_gatk=False``.
        hq_sites_gq_thresh (int): For finding high quality sites, the minimum GQ threshold for a kept genotype.
        hq_sites_ab_thresh (float): For finding high quality sites, the allowable allelic balance range for a kept heterozygote genotype.
        hq_sites_maf_thresh (float): For finding high quality sites, the minimum minor allele frequency as computed by `hl.variant_qc` in order for a variant to be considered to be common.
        hq_sites_mac_thresh (int): For finding high quality sites, the minimum number of minor alleles observed as computed by `hl.variant_qc` in order for a variant to be considered common.
        hq_sites_call_rate_thresh (float): For finding high quality sites, the minimum call rate of a variant as computed by `hl.variant_qc` in order for a variant to be considered.
        pre_pruned_variants (Optional(Table)): For subsetting data to LD pruned sites, the list of pre-pruned variants to use.
        r2_thresh (float): For subsetting data to LD pruned sites, the R2 threshold to use when running ``hl.ld_prune``.
        bp_window_size (int): For subsetting data to LD pruned sites, the window for considering variants to prune in base pairs.
        sex_reported_is_female (Optional(BooleanExpression)): For running a sex check, the reported sex as a Boolean Expression keyed by the Sample ID where ``True`` is Female and ``False`` is Male.
        sex_check_male_fhet_thresh (float): For running a sex check, the minimum Fstat threshold for calling a sample "Male".
        sex_check_female_fhet_thresh (float): For running a sex check, the maximum Fstat threshold for calling a sample "Female".
        pcs_k (int): For computing principal components, the number of principal components to compute.
        ancestry_pca_loadings (str): For computing ancestry pop labels, path to the gnomAD PCA loadings Hail Table.
        ancestry_onnx_rf (str): For computing ancestry pop labels, path to the gnomAD random forest model.
        ancestry_num_pcs (int): For computing ancestry pop labels, number of PCs to use when applying the random forest model. This is dependent on the onnx model.
        ancestry_min_prob (float): For computing ancestry pop labels, the minimum probability for a predicted population before giving a sample that label.
        chimera_rate (Optional(NumericExpression)): For filtering samples based on chimera rate, the chimera rate as a NumericExpression keyed by the sample ID.
        chimera_threshold (float): For filtering samples based on chimera rate, the maximum threshold to use at which a sample is considered passing.
        contamination_charr_thresh (float): For filtering samples based on contamination rate, the threshold for which a sample is considered to be passing from CHARR.
        contamination_min_af (float): For filtering samples based on contamination rate, the minimum allele frequency when using the ``compute_charr`` function in Hail.
        contamination_max_af (float): For filtering samples based on contamination rate, the maximum allele frequency when using the ``compute_charr`` function in Hail.
        contamination_min_dp (int): For filtering samples based on contamination rate, the minimum depth when using the ``compute_charr`` function in Hail.
        contamination_max_dp (int): For filtering samples based on contamination rate, the maximum depth when using the ``compute_charr`` function in Hail.
        contamination_min_gq (int): For filtering samples based on contamination rate, the minimum GQ when using the ``compute_charr`` function in Hail.
        contamination_ref_af (Optional[NumericExpression]): A float row field on the MatrixTable with the reference allele frequency when using the ``compute_charr`` function in Hail.
        coverage_gq_thresh (int): For filtering samples based on GQ coverage stats, the minimum GQ for a genotype to be considered passing.
        coverage_fraction_thresh (float): For filtering samples based on GQ coverage stats, the minimum rate at which GQ is above the minimum threshold.
        log_path (Optional[str]): Log file path to write to.
        hail_init_kwargs (Optional[dict]): Keyword arguments to hl.init().
        log (Optional[logging.Logger]): Logging object to log to.
        make_plots (bool): Save default QC plots to ``out_dir``.

    Returns:
        A :class:`.SampleQCResult` object.
    """
    out_dir = out_dir.rstrip('/') + '/'

    if log is None:
        log = setup_logger()

    if isinstance(dataset, hl.MatrixTable):
        init_n_variants, init_n_samples = mt.count()
    else:
        assert isinstance(dataset, hl.vds.VariantDataset)
        init_n_variants, init_n_samples = dataset.variant_data.count()

    log.info(f'QC of dataset with {init_n_variants} variants and {init_n_samples} samples.')

    if passing_variants is None:
        log.info(f'selecting passing variants based on filters')
        passing_variants = select_passing_variants_from_filters(dataset)

    log.info(f'''selecting high quality common sites (`select_high_quality_common_sites`) with parameters:
is_gatk={is_gatk}
dp_thresh={hq_sites_dp_thresh}
gq_thresh={hq_sites_gq_thresh}
ab_thresh={hq_sites_ab_thresh}
maf_thresh={hq_sites_maf_thresh}
mac_thresh={hq_sites_mac_thresh}
call_rate_thresh={hq_sites_call_rate_thresh}
''')

    high_qual_sites = select_high_quality_common_sites(dataset=dataset,
                                                       is_gatk=is_gatk,
                                                       exome_regions=exome_regions,
                                                       low_complexity_regions=low_complexity_regions,
                                                       passing_variants=passing_variants,
                                                       dp_thresh=hq_sites_dp_thresh,
                                                       gq_thresh=hq_sites_gq_thresh,
                                                       ab_thresh=hq_sites_ab_thresh,
                                                       maf_thresh=hq_sites_maf_thresh,
                                                       mac_thresh=hq_sites_mac_thresh,
                                                       call_rate_thresh=hq_sites_call_rate_thresh)

    high_qual_sites = high_qual_sites.persist()

    log.info(f'''inferring ancestry population labels (`infer_ancestry`) with parameters:
is_gatk={is_gatk}
dp_thresh={hq_sites_dp_thresh}
gq_thresh={hq_sites_gq_thresh}
ab_thresh={hq_sites_ab_thresh}
pca_loadings={ancestry_pca_loadings}
onnx_rf={ancestry_onnx_rf}
num_pcs={ancestry_num_pcs}
min_prob={ancestry_min_prob}
''')

    ancestry_pop_labels = infer_ancestry(dataset=dataset,
                                         is_gatk=is_gatk,
                                         dp_thresh=hq_sites_dp_thresh,
                                         gq_thresh=hq_sites_gq_thresh,
                                         ab_thresh=hq_sites_ab_thresh,
                                         pca_loadings=ancestry_pca_loadings,
                                         onnx_rf=ancestry_onnx_rf,
                                         num_pcs=ancestry_num_pcs,
                                         min_prob=ancestry_min_prob,
                                         log=log)

    log.info(f'''computing sample qc statistics (`calculate_sample_qc_stats`) with parameters:
sex_female_fhet_thresh={sex_check_female_fhet_thresh}
sex_male_fhet_thresh={sex_check_male_fhet_thresh}
pcs_k={pcs_k}
chimera_threshold={chimera_threshold},
contamination_charr_thresh={contamination_charr_thresh}
contamination_min_af={contamination_min_af}
contamination_max_af={contamination_max_af}
contamination_min_dp={contamination_min_dp}
contamination_max_dp={contamination_max_dp}
contamination_min_gq={contamination_min_gq}
coverage_gq_thresh={coverage_gq_thresh}
coverage_fraction_thresh={coverage_fraction_thresh}
ld_prune_r2_thresh={r2_thresh}
ld_prune_bp_window_size={bp_window_size}
''')

    sample_qc_stats = calculate_sample_qc_stats(dataset=dataset,
                                                high_quality_sites=high_qual_sites,
                                                sex_female_fhet_thresh=sex_check_female_fhet_thresh,
                                                sex_male_fhet_thresh=sex_check_male_fhet_thresh,
                                                sex_reported_is_female=sex_reported_is_female,
                                                # relatedness_min_kinship=relatedness_min_kinship,
                                                pcs_k=pcs_k,
                                                chimera_rate=chimera_rate,
                                                chimera_threshold=chimera_threshold,
                                                contamination_charr_thresh=contamination_charr_thresh,
                                                contamination_min_af=contamination_min_af,
                                                contamination_max_af=contamination_max_af,
                                                contamination_min_dp=contamination_min_dp,
                                                contamination_max_dp=contamination_max_dp,
                                                contamination_min_gq=contamination_min_gq,
                                                contamination_ref_AF=contamination_ref_af,
                                                coverage_gq_thresh=coverage_gq_thresh,
                                                coverage_fraction_thresh=coverage_fraction_thresh,
                                                ld_prune_r2_thresh=r2_thresh,
                                                ld_prune_bp_window_size=bp_window_size,
                                                pre_pruned_variants=pre_pruned_variants,
                                                ancestry_pop_labels=ancestry_pop_labels,
                                                log=log)

    sample_qc_stats.write(f'{out_dir}sample_qc_stats.ht', overwrite=True)
    log.info(f'wrote out raw sample QC statistics to {out_dir}sample_qc_stats.ht')

    sample_qc_stats = SampleQCResult.load(f'{out_dir}sample_qc_stats.ht')

    if make_plots:
        log.info(f'writing out default plots to {out_dir}')

        pc1_pc2 = sample_qc_stats.plot_pcs(1, 2)
        save_ggplot(pc1_pc2, f'{out_dir}pcs/pc1_pc2.png')

        pc1_pc3 = sample_qc_stats.plot_pcs(1, 3)
        save_ggplot(pc1_pc3, f'{out_dir}pcs/pc1_pc3.png')

        pc2_pc3 = sample_qc_stats.plot_pcs(2, 3)
        save_ggplot(pc2_pc3, f'{out_dir}pcs/pc2_pc3.png')

        box_plots = sample_qc_stats.plot_qc_metric_boxplots()
        for metric, plot in box_plots.items():
            save_ggplot(plot, f'{out_dir}qc/{metric}_boxplot.png')

        density_plots = sample_qc_stats.plot_qc_metric_density()
        for metric, plot in density_plots.items():
            save_ggplot(plot, f'{out_dir}qc/{metric}_density.png')

        passing_samples = sample_qc_stats.passing_samples()
        passing_samples.sample_id.export(f'{out_dir}passing_sample_ids.tsv', header=True)

        pass_box_plots = passing_samples.plot_qc_metric_boxplots()
        for metric, plot in pass_box_plots.items():
            save_ggplot(plot, f'{out_dir}qc/{metric}_pass_boxplot.png')

        pass_density_plots = passing_samples.plot_qc_metric_density()
        for metric, plot in pass_density_plots.items():
            save_ggplot(plot, f'{out_dir}qc/{metric}_pass_density.png')

    return sample_qc_stats
