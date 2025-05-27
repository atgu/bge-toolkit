import logging
from typing import Optional, Tuple

import hail as hl
from hail.backend.service_backend import ServiceBackend

from ..common.io import load_input_data
from ..common.binned_exprs import (mac_bin, maf_bin, max_gp_bin, gq_bin, qual_approx_bin, info_score_bin, dp_bin)
from ..common.utils import setup_logger
from .joint_callset import JointCallSet
from .concordance import Statistic
from .sample_qc import calculate_sample_qc_stats, select_high_quality_sites
from .utils import ALL_AGG, ALL_GROUP, JoinType


def setup_logging_and_qob(log_name: str,
                          log_path: Optional[str] = None,
                          hail_init_kwargs: Optional[dict] = None,
                          log: Optional[logging.Logger] = None) -> logging.Logger:
    log = setup_logger(log, log_path, log_name)
    if log is None:
        log = logging.getLogger(log_name)
        log.propagate = False

    bad_logs = ['hailtop.aiocloud.aiogoogle.credentials',
                'batch_client.aioclient']

    for bad_log in bad_logs:
        logging.getLogger(bad_log).setLevel(logging.WARNING)

    if hail_init_kwargs:
        hl.init(**hail_init_kwargs)

    backend = hl.current_backend()
    assert isinstance(backend, ServiceBackend)
    backend.disable_progress_bar = False

    return log


def load_data(*, path: str, descriptor: str, log: logging.Logger, n_partitions: Optional[int] = None) -> hl.MatrixTable:
    log.info(f'loading dataset {descriptor}')
    mt = load_input_data(path, n_partitions=n_partitions)
    n_variants, n_samples = mt.count()
    log.info(f'found {n_variants} variants and {n_samples} samples in {descriptor}.')
    return mt


def apply_filters(*,
                  mt: hl.MatrixTable,
                  description: str,
                  log: logging.Logger,
                  sample_list: Optional[str] = None,
                  variant_list: Optional[str] = None,
                  contig: Optional[str] = None,
                  n_samples: Optional[int] = None,
                  n_variants: Optional[int] = None,
                  downsample_samples: Optional[float] = None,
                  downsample_variants: Optional[float] = None) -> Tuple[hl.MatrixTable, bool]:
    downsampled = False

    if contig is not None:
        downsampled = True
        mt = mt.filter_rows(mt.locus.contig == contig)
        log.info(f'Subsetting to variants on {contig} for {description}.')

    if downsample_variants:
        downsampled = True
        mt = mt.sample_rows(downsample_variants, seed=10243523)
        log.info(f'Downsampling variants by {downsample_variants} for {description}.')

    if downsample_samples:
        downsampled = True
        mt = mt.sample_cols(downsample_samples, seed=10243523)
        log.info(f'Downsampling samples by {downsample_samples} for {description}.')

    if n_variants is not None:
        downsampled = True
        mt = mt.head(n_rows=n_variants)
        log.info(f'Selecting the first {n_variants} variants for {description}.')

    if n_samples is not None:
        downsampled = True
        mt = mt.head(n_rows=None, n_cols=n_samples)
        log.info(f'Selecting the first {n_samples} samples for {description}.')

    if sample_list is not None:
        samples = hl.import_table(sample_list)
        samples = samples.key_by('s')
        mt = mt.semi_join_cols(samples)
        log.info(f'Subset to samples in {sample_list} for {description}.')

    if variant_list is not None:
        downsampled = True
        variants = hl.import_table(variant_list, no_header=True)
        variants = variants.key_by('v')
        variants = variants.annotate(variant=hl.parse_variant(variants.v))
        variants = variants.annotate(locus=variants.variant.locus, alleles=variants.variant.alleles)
        variants = variants.key_by('locus', 'alleles')

        mt = mt.semi_join_rows(variants)

        log.info(f'Subset to variants in {variant_list} for {description}.')

    return (mt, downsampled)


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
                contig: Optional[str] = None,
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
        contig (str): Filter to variants in the contig.
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


def sample_qc_gatk(path: str,
                   out_dir: str,
                   exome_regions: Optional[str],
                   low_complexity_regions: Optional[str],
                   chimeric_read_rate: Optional[Tuple[str, str, str, float]],
                   relatedness_pi_hat: float,
                   outlier_z_score: int,
                   maf: float,
                   mac: int,
                   prune_r2: float,
                   read_depth: Tuple[float, float],
                   male_f_het: float,
                   female_f_het: float,
                   sample_list: Optional[str],
                   variant_list: Optional[str],
                   contig: Optional[str],
                   n_samples: Optional[int],
                   n_variants: Optional[int],
                   downsample_samples: Optional[float],
                   downsample_variants: Optional[float],
                   log_path: Optional[str],
                   hail_init_kwargs: Optional[dict],
                   log: Optional[logging.Logger],
                   n_partitions: Optional[int],
                   reference_genome: Optional[str]):
    log = setup_logging_and_qob(log_name='sample_qc_gatk',
                                log_path=log_path,
                                hail_init_kwargs=hail_init_kwargs,
                                log=log)

    mt = load_data(path=path, descriptor=path, log=log, n_partitions=n_partitions)

    mt, _ = apply_filters(mt=mt,
                          description=path,
                          log=log,
                          sample_list=sample_list,
                          variant_list=variant_list,
                          contig=contig,
                          n_samples=n_samples,
                          n_variants=n_variants,
                          downsample_samples=downsample_samples,
                          downsample_variants=downsample_variants)

    out_dir = out_dir.rstrip('/') + '/'

    exome_regions = hl.import_locus_intervals(exome_regions, reference_genome=reference_genome)
    low_complexity_regions = hl.import_locus_intervals(low_complexity_regions, reference_genome=reference_genome)

    high_qual_sites = select_high_quality_sites(mt,
                                                exome_regions=exome_regions,
                                                low_complexity_regions=low_complexity_regions)

    sample_qc_stats = calculate_sample_qc_stats(mt, high_quality_sites=high_qual_sites)

    # Ultimate QC filter condition based on all of the stats computed
    passes_qc = ...

    sample_qc_stats = sample_qc_stats.annotate(passes_qc=passes_qc)

    sample_qc_stats.write(f'{out_dir}/sample_qc_stats.ht', overwrite=True)

    passing_samples = sample_qc_stats.filter(sample_qc_stats.passes_qc).select()
    passing_samples.export(f'{out_dir}/passing_samples.tsv', sep="\t")
