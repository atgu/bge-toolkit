from dataclasses import asdict
from typing import Annotated as Ann, Tuple

import typer
from typer import Option as Opt

from .utils import JoinType
from .main import concordance as _concordance, sample_qc_gatk as _sample_qc_gatk


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

sample_qc_gatk_app = typer.Typer(
    name='sample-qc-gatk',
    invoke_without_command=True,
    no_args_is_help=True,
    help='Retain high quality samples from a GATK called exome dataset.',
)

qc_app.add_typer(concordance_app)
qc_app.add_typer(sample_qc_gatk_app)


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
                contig: Ann[str, Opt('--contig', help='Filter to variants in the contig.')] = None,
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


@sample_qc_gatk_app.callback()
def sample_qc_gatk(ctx: typer.Context,
                   exome_path: Ann[str, Opt('--exome', help='Exome dataset called with GATK to filter out low quality samples.')],
                   out_dir: Ann[str, Opt('--output-dir', help='Output directory.')],
                   exome_regions: Ann[str, Opt('--exome-regions', help='BED file with list of exome intervals.')] = None,
                   chimeric_read_rate: Ann[Tuple[str, str, str, float], Opt('--chimeric-read-rate', help='File with chimeric read rates along with the sample ID column name, the chimeric read rate column name, and the threshold for filtering')] = None,
                   relatedness_pi_hat: Ann[float, Opt('--relatedness-pi-hat', help='Threshold for related individuals.')] = 0.1,
                   outlier_z_score: Ann[int, Opt('--outlier-z-score', help='z-score for considering a sample an outlier on QC metrics.')] = 5,
                   maf: Ann[float, Opt('--maf', help='MAF threshold for common variants.')] = 0.01,
                   mac: Ann[int, Opt('--mac', help='Minimum MAC threshold')] = 10,
                   prune_r2: Ann[float, Opt('--prune-r2', help='LD-prune R2 value for pruning common variants.')] = 0.2,
                   read_depth: Ann[Tuple[float, float], Opt('--read-depth', help='Two values representing the percentage of the genome (first value) with read depth at least a given threshold (second value).')] = (0.9, 20),
                   male_f_het: Ann[float, Opt('--male-f-het', help='Fhet threshold for calling a sample a male.')] = 0.8,
                   female_f_het: Ann[float, Opt('--female-f-het', help='Fhet threshold for calling a sample a female.')] = 0.2,
                   sample_list: Ann[str, Opt('--sample-list', help='Filter to samples listed in the file.')] = None,
                   variant_list: Ann[str, Opt('--variant-list', help='Filter to variants listed in the file.')] = None,
                   contig: Ann[str, Opt('--contig', help='Filter to variants in the contig.')] = None,
                   n_samples: Ann[int, Opt('--n-samples', help='Filter to the first N overlapping samples.')] = None,
                   n_variants: Ann[int, Opt('--n-variants', help='Filter to the first N variants.')] = None,
                   downsample_samples: Ann[float, Opt('--downsample-samples', help='Downsample to X fraction of samples.')] = None,
                   downsample_variants: Ann[float, Opt('--downsample-variants', help='Downsample to X fraction of variants.')] = None,
                   log_path: Ann[str, Opt('--log-path', help='Log file path to write to.')] = None,
                   n_exome_partitions: Ann[int, Opt('--n-exome-partitions', help="Number of partitions to load the exome dataset with if it's a MatrixTable.")] = None,
                   ):
    _sample_qc_gatk(exome_path,
                    out_dir,
                    exome_regions,
                    chimeric_read_rate,
                    relatedness_pi_hat,
                    ...)
