from dataclasses import asdict
from typing import Annotated as Ann, List

import typer
from typer import Option as Opt

import hail as hl

from .utils import JoinType
from .main import concordance as _concordance, sample_qc as _sample_qc
from ..common.io import load_data
from ..common.utils import apply_filters, setup_logging_and_qob

qc_app = typer.Typer(
    name='qc',
    no_args_is_help=True,
    help='Run QC-related tools on BGE datasets.',
    pretty_exceptions_show_locals=False,
)

concordance_app = typer.Typer(
    name='concordance',
    invoke_without_command=True,
    no_args_is_help=True,
    help='Run concordance on BGE exome and imputation datasets.',
)

sample_qc_app = typer.Typer(
    name='sample-qc',
    invoke_without_command=True,
    no_args_is_help=True,
    help='Compute sample QC statistics from an exome dataset.',
)

qc_app.add_typer(concordance_app)
qc_app.add_typer(sample_qc_app)


@concordance_app.callback()
def concordance(ctx: typer.Context,
                exome_path: Ann[str, Opt('--exome', help='Exome dataset to compare to.')],
                imputation_path: Ann[str, Opt('--imputation', help='Imputation dataset to compare to.')],
                out_dir: Ann[str, Opt('--output-dir', help='Output directory.')],
                exome_maf: Ann[bool, Opt('--EXOME-MAF', help='Bin concordance counts by Exome Minor Allele Frequency. (Global+Variant)')] = False,
                exome_mac: Ann[bool, Opt('--EXOME-MAC', help='Bin concordance counts by Exome Minor Allele Counts. (Global+Variant)')] = False,
                imputation_maf: Ann[bool, Opt('--IMPUTATION-MAF', help='Bin concordance counts by Imputation Minor Allele Frequency. (Global+Variant)')] = False,
                imputation_mac: Ann[bool, Opt('--IMPUTATION-MAC', help='Bin concordance counts by Imputation Minor Allele Counts. (Global+Variant)')] = False,
                info_score: Ann[bool, Opt('--INFO', help='Bin concordance counts by imputation INFO score. (Global+Variant)')] = False,
                dp: Ann[bool, Opt('--DP', help='Bin concordance counts by exome genotype DP. (Global+Variant)')] = False,
                gq: Ann[bool, Opt('--GQ', help='Bin concordance counts by exome genotype GQ. (Global+Variant)')] = False,
                max_gp: Ann[bool, Opt('--MAX-GP', help='Bin concordance counts by maximum value of imputation GP. (Global+Variant)')] = False,
                qual_approx: Ann[bool, Opt('--QUAL-APPROX', help='Bin concordance counts by variant Approx Qual Score. (Global+Variant)')] = False,
                sample_list: Ann[str, Opt('--sample-list', help='Filter to samples listed in the file.')] = None,
                variant_list: Ann[str, Opt('--variant-list', help='Filter to variants listed in the file.')] = None,
                contig: Ann[List[str], Opt('--contig', help='Filter to variants in the contig.')] = None,
                n_samples: Ann[int, Opt('--n-samples', help='Filter to the first N overlapping samples.')] = None,
                n_variants: Ann[int, Opt('--n-variants', help='Filter to the first N variants.')] = None,
                downsample_samples: Ann[float, Opt('--downsample-samples', help='Downsample to X fraction of samples.')] = None,
                downsample_variants: Ann[float, Opt('--downsample-variants', help='Downsample to X fraction of variants.')] = None,
                join_type: Ann[JoinType, Opt('--join-type', show_choices=True, help='Join type.')] = JoinType.inner,
                run_variants: Ann[bool, Opt('--variant-conc', help='Run variant concordance.')] = False,
                run_samples: Ann[bool, Opt('--sample-conc', help='Run sample concordance.')] = False,
                run_global: Ann[bool, Opt('--global-conc', help='Run global concordance.')] = False,
                log_path: Ann[str, Opt('--log-path', help='Log file path to write to.')] = None,
                n_exome_partitions: Ann[int, Opt('--n-exome-partitions', help="Number of partitions to load the exome dataset with if it's a MatrixTable.")] = None,
                n_imp_partitions: Ann[int, Opt('--n-imp-partitions', help="Number of partitions to load the imputation dataset with if it's a MatrixTable.")] = None,
                ):
    _concordance(exome_path=exome_path,
                 imputation_path=imputation_path,
                 out_dir=out_dir,
                 exome_maf=exome_maf,
                 exome_mac=exome_mac,
                 imputation_maf=imputation_maf,
                 imputation_mac=imputation_mac,
                 info_score=info_score,
                 dp=dp,
                 gq=gq,
                 max_gp=max_gp,
                 qual_approx=qual_approx,
                 sample_list=sample_list,
                 variant_list=variant_list,
                 contig=contig,
                 n_samples=n_samples,
                 n_variants=n_variants,
                 downsample_samples=downsample_samples,
                 downsample_variants=downsample_variants,
                 join_type=join_type,
                 run_variants=run_variants,
                 run_samples=run_samples,
                 run_global=run_global,
                 log_path=log_path,
                 hail_init_kwargs=asdict(ctx.obj),
                 log=None,
                 n_exome_partitions=n_exome_partitions,
                 n_imp_partitions=n_imp_partitions)


@sample_qc_app.callback()
def sample_qc(ctx: typer.Context,
              dataset_path: Ann[str, Opt('--dataset', help='Path to MatrixTable or VariantDataset.')],
              out_dir: Ann[str, Opt('--output-dir', help='Output directory.')],
              exome_regions: Ann[str, Opt('--exome-regions', help='BED file with list of exome intervals.')],
              low_complexity_regions: Ann[str, Opt('--low-complexity-regions', help='BED file with list of exome intervals.')],
              passing_variants: Ann[str, Opt('--passing-variants', help='Path to Hail Table with list of passing variants from VQSR if available')] = None,
              gatk: Ann[bool, Opt('--gatk', help='GATK was used to generate the callset.')] = False,
              dragen: Ann[bool, Opt('--dragen', help='DRAGEN was used to generate the callset.')] = False,
              hq_sites_dp_thresh: Ann[int, Opt('--hq-sites-dp-thresh', help='DP threshold for marking a genotype high quality.')] = 10,
              hq_sites_gq_thresh: Ann[int, Opt('--hq-sites-gq-thresh', help='GQ threshold for marking a genotype high quality.')] = 20,
              hq_sites_ab_thresh: Ann[int, Opt('--hq-sites-ab-thresh', help='AB threshold for marking a heterozygous genotype high quality.')] = 0.2,
              hq_sites_maf_thresh: Ann[int, Opt('--hq-sites-maf-thresh', help='MAF threshold for marking a variant high quality.')] = 0.01,
              hq_sites_mac_thresh: Ann[int, Opt('--hq-sites-mac-thresh', help='MAC threshold for marking a variant high quality.')] = 10,
              hq_sites_call_rate_thresh: Ann[int, Opt('--hq-sites-call-rate-thresh', help='Call rate threshold for marking a variant high quality.')] = 0.95,
              pre_pruned_variants: Ann[str, Opt('--pre-pruned-variants', help='Path to Hail Table with variant sites in the same reference genome as the exome data.')] = None,
              r2_thresh: Ann[float, Opt('--r2-thresh', help='r2 threshold for pruning variants when computing PCs and relatedness.')] = 0.2,
              bp_window_size: Ann[int, Opt('--bp-window-size', help='BP window size for pruning variants when computing PCs and relatedness')] = 1000000,
              reported_sex_table_path: Ann[str, Opt('--reported-sex-path', help='Path to Hail table where reported sex field is located.')] = None,
              reported_sex_col_name: Ann[str, Opt('--reported-sex-col', help='Column name with the reported sex coded as True for female and False for male')] = None,
              sex_check_male_fhet_thresh: Ann[float, Opt('--sex-check-male-fhet-thresh', help='Sex check male Fhet threshold.')] = 0.8,
              sex_check_female_fhet_thresh: Ann[float, Opt('--sex-check-female-fhet-thresh', help='Sex check female Fhet threshold.')] = 0.2,
              # relatedness_min_kinship: Ann[float, Opt('--relatedness-min-kinship', help='Relatedness min kinship to be considered related.')] = 0.125,
              pcs_k: Ann[int, Opt('--pcs-k', help='Number of principal components to compute.')] = 10,
              ancestry_pca_loadings: Ann[str, Opt('--ancestry-pop-loadings', help='path to gnomAD loadings table')] = 'gs://gcp-public-data--gnomad/release/4.0/pca/gnomad.v4.0.pca_loadings.ht',
              ancestry_onnx_rf: Ann[str, Opt('--ancestry-onnx-rf', help='path to gnomAD onnx file')] = 'gs://gcp-public-data--gnomad/release/4.0/pca/gnomad.v4.0.RF_fit.onnx',
              ancestry_num_pcs: Ann[int, Opt('--ancestry-pcs-k', help='Number of principal components to use. This is dependent on the gnomAD data.')] = 20,
              ancestry_min_prob: Ann[float, Opt('--ancestry-min-prob', help='Minimum probability to be called for an ancestry population.')] = 0.75,
              chimera_rate_path: Ann[str, Opt('--chimera-rate-path', help='Path to Hail Table with chimera read rates.')] = None,
              chimera_rate_col: Ann[str, Opt('--chimera-rate-col', help='Column name of chimera read rates in the Hail Table.')] = None,
              chimera_threshold: Ann[float, Opt('--chimera-threshold', help='Maximum rate of chimera reads for a sample to be considered passing.')] = 0.05,
              contamination_charr_thresh: Ann[float, Opt('--contamination-charr-thresh', help='Contamination rate threshold as computed using CHARR.')] = 0.05,
              contamination_min_af: Ann[float, Opt('--contamination-min-af', help='Minimum AF when computing contamination rating using CHARR.')] = 0.05,
              contamination_max_af: Ann[float, Opt('--contamination-max-af', help='Maximum AF when computing contamination rating using CHARR.')] = 0.95,
              contamination_min_dp: Ann[int, Opt('--contamination-min-dp', help='Minimum DP when computing contamination rating using CHARR.')] = 10,
              contamination_max_dp: Ann[int, Opt('--contamination-max-dp', help='Maximum DP when computing contamination rating using CHARR.')] = 100,
              contamination_min_gq: Ann[int, Opt('--contamination-min-gq', help='Minimum GQ when computing contamination rating using CHARR.')] = 20,
              contamination_ref_af_path: Ann[str, Opt('--contamination-ref-af-path', help='Path to Hail Table with contamination rates per sample.')] = None,
              contamination_ref_af_col_name: Ann[str, Opt('--contamination-ref-af-col-name', help='Column name of contamination rate in Hail Table.')] = None,
              coverage_gq_thresh: Ann[int, Opt('--coverage-gq-thresh', help='Minimum GQ for considering a variant well-genotyped in an individual.')] = 20,
              coverage_fraction_thresh: Ann[float, Opt('--coverage-fraction-thresh', help='Minimum number of genotypes meeting the GQ threshold for being considered passing.')] = 0.9,
              sample_list: Ann[str, Opt('--sample-list', help='Filter to samples listed in the file.')] = None,
              variant_list: Ann[str, Opt('--variant-list', help='Filter to variants listed in the file.')] = None,
              contig: Ann[List[str], Opt('--contig', help='Filter to variants in the contig.')] = None,
              n_samples: Ann[int, Opt('--n-samples', help='Filter to the first N overlapping samples.')] = None,
              n_variants: Ann[int, Opt('--n-variants', help='Filter to the first N variants.')] = None,
              downsample_samples: Ann[float, Opt('--downsample-samples', help='Downsample to X fraction of samples.')] = None,
              downsample_variants: Ann[float, Opt('--downsample-variants', help='Downsample to X fraction of variants.')] = None,
              log_path: Ann[str, Opt('--log-path', help='Log file path to write to.')] = None,
              n_partitions: Ann[int, Opt('--n-partitions', help="Number of partitions to load the dataset with if it's a MatrixTable.")] = None,
              reference_genome: Ann[str, Opt('--reference-genome', help='Reference genome string to use when loading BED file intervals.')] = None):
    log = setup_logging_and_qob(log_name='sample_qc',
                                log_path=log_path,
                                hail_init_kwargs=asdict(ctx.obj),
                                log=None)

    if gatk and dragen:
        raise ValueError('cannot specify both --gatk and --dragen')
    if not gatk and not dragen:
        raise ValueError('must specify one of --gatk or --dragen')

    if pre_pruned_variants is not None:
        pre_pruned_variants = hl.read_table(pre_pruned_variants)

    if reported_sex_table_path is not None:
        if reported_sex_col_name is None:
            raise ValueError(f'must specify the column name with the boolean reported sex field in {reported_sex_table_path}')
        reported_sex = hl.read_table(reported_sex_table_path)[reported_sex_col_name]
    else:
        reported_sex = None

    if chimera_rate_path is not None:
        if chimera_rate_col is None:
            raise ValueError(f'must specify the column name with the chimera rate field in {chimera_rate_path}')
        chimera_rate = hl.read_table(chimera_rate_path)[chimera_rate_col]
    else:
        chimera_rate = None

    if contamination_ref_af_path is not None:
        if contamination_ref_af_col_name is None:
            raise ValueError(f'must specify the column name with the contamination rate field in {contamination_ref_af_path}')
        contamination_ref_af = hl.read_table(contamination_ref_af_path)[contamination_ref_af_col_name]
    else:
        contamination_ref_af = None

    if passing_variants is not None:
        passing_variants = hl.read_table(passing_variants)

    dataset = load_data(path=dataset_path, descriptor='dataset', log=log, n_partitions=n_partitions)

    dataset, _ = apply_filters(dataset=dataset_path,
                               description=dataset_path,
                               log=log,
                               sample_list=sample_list,
                               variant_list=variant_list,
                               contig=contig,
                               n_samples=n_samples,
                               n_variants=n_variants,
                               downsample_samples=downsample_samples,
                               downsample_variants=downsample_variants)

    if reference_genome is None:
        reference_genome = dataset.locus.dtype.reference_genome

    exome_regions = hl.import_locus_intervals(exome_regions, reference_genome=reference_genome)
    low_complexity_regions = hl.import_locus_intervals(low_complexity_regions, reference_genome=reference_genome, skip_invalid_intervals=True)

    _sample_qc(dataset=dataset,
               out_dir=out_dir,
               is_gatk=gatk,
               exome_regions=exome_regions,
               low_complexity_regions=low_complexity_regions,
               passing_variants=passing_variants,
               hq_sites_dp_thresh=hq_sites_dp_thresh,
               hq_sites_gq_thresh=hq_sites_gq_thresh,
               hq_sites_ab_thresh=hq_sites_ab_thresh,
               hq_sites_maf_thresh=hq_sites_maf_thresh,
               hq_sites_mac_thresh=hq_sites_mac_thresh,
               hq_sites_call_rate_thresh=hq_sites_call_rate_thresh,
               pre_pruned_variants=pre_pruned_variants,
               r2_thresh=r2_thresh,
               bp_window_size=bp_window_size,
               sex_reported_is_female=reported_sex,
               sex_check_male_fhet_thresh=sex_check_male_fhet_thresh,
               sex_check_female_fhet_thresh=sex_check_female_fhet_thresh,
               # relatedness_min_kinship=relatedness_min_kinship,
               pcs_k=pcs_k,
               ancestry_pca_loadings=ancestry_pca_loadings,
               ancestry_onnx_rf=ancestry_onnx_rf,
               ancestry_num_pcs=ancestry_num_pcs,
               ancestry_min_prob=ancestry_min_prob,
               chimera_rate=chimera_rate,
               chimera_threshold=chimera_threshold,
               contamination_charr_thresh=contamination_charr_thresh,
               contamination_min_af=contamination_min_af,
               contamination_max_af=contamination_max_af,
               contamination_min_dp=contamination_min_dp,
               contamination_max_dp=contamination_max_dp,
               contamination_min_gq=contamination_min_gq,
               contamination_ref_af=contamination_ref_af,
               coverage_gq_thresh=coverage_gq_thresh,
               coverage_fraction_thresh=coverage_fraction_thresh,
               log=log,
               make_plots=True)
