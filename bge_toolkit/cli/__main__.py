import json
from dataclasses import dataclass
from typing import Annotated as Ann, Dict, List, Optional, Tuple, Union

import typer
from typer import Option as Opt

from ..qc import cli as qc_cli


app = typer.Typer(
    help='Command line interface for the BGE-Toolkit.',
    no_args_is_help=True,
    pretty_exceptions_show_locals=False,
)


@dataclass
class HailInitKwargs:
    app_name: Optional[str]
    master: Optional[str]
    local: str
    log: Optional[str]
    quiet: bool
    append: bool
    min_block_size: int
    branching_factor: int
    tmp_dir: Optional[str]
    default_reference: Optional[str]
    idempotent: bool
    global_seed: Optional[int]
    spark_conf: Optional[Dict[str, str]]
    skip_logging_configuration: bool
    local_tmpdir: Optional[str]
    _optimizer_iterations: Optional[int]
    backend: Optional[str]
    driver_cores: Optional[int]
    driver_memory: Optional[str]
    worker_cores: Optional[int]
    worker_memory: Optional[str]
    gcs_requester_pays_configuration: Optional[Union[str, Tuple[str, List[str]]]]
    regions: Optional[List[str]]
    gcs_bucket_allow_list: Optional[Dict[str, List[str]]]
    copy_spark_log_on_error: Optional[bool]


@app.callback()
def common(ctx: typer.Context,
           app_name: Ann[str, Opt('--app-name', help='A name for this pipeline. In the Spark backend, this becomes the Spark application name. '
                                             'In the Batch backend, this is a prefix for the name of every Batch.')] = None,
           master: Ann[str, Opt('--master', help='Spark Backend only. URL identifying the Spark leader (master) node or `local[N]` for local clusters.')] = None,
           local: Ann[str, Opt('--local', help='Spark Backend only. Local-mode core limit indicator. Must either be `local[N]` where N is a '
                                               'positive integer or `local[*]`. The latter indicates Spark should use all cores '
                                               'available. `local[*]` does not respect most containerization CPU limits. This option is only '
                                               'used if `master` is unset and `spark.master` is not set in the Spark configuration.')] = 'local[*]',
           log: Ann[str, Opt('--log', help='Local path for Hail log file. Does not currently support distributed file systems like Google Storage, S3, or HDFS.')] = None,
           quiet: Ann[bool, Opt('--quiet', help='Print fewer log messages.')] = False,
           append: Ann[bool, Opt('--append', help='Append to the end of the log file.')] = False,
           min_block_size: Ann[int, Opt('--min-block-size', help='Minimum file block size in MB.')] = 0,
           branching_factor: Ann[int, Opt('--branching-factor', help='Branching factor for tree aggregation.')] = 50,
           tmp_dir: Ann[str, Opt('--tmp-dir', help='Networked temporary directory. Must be a network-visible file path. Defaults to /tmp in the default scheme.')] = None,
           default_reference: Ann[str, Opt('--default-reference', help='Default reference genome.')] = None,
           idempotent: Ann[bool, Opt('--idempotent', help='If ``True``, calling this function is a no-op if Hail has already been initialized.')] = False,
           global_seed: Ann[int, Opt('--global-seed', help='Global random seed.')] = None,
           spark_conf: Ann[str, Opt('--spark-conf', help='Spark backend only. Spark configuration parameters in json format.')] = None,
           skip_logging_configuration: Ann[bool, Opt('--skip-logging-configuration', help='Spark Backend only. Skip logging configuration in java and python.')] = False,
           local_tmpdir: Ann[str, Opt('--local-tmpdir', help='Local temporary directory. Used on driver and executor nodes. Must use the file scheme. Defaults to TMPDIR, or /tmp.')] = None,
           _optimizer_iterations: Ann[int, Opt('--optimizer-iterations', help='')] = None,
           backend: Ann[str, Opt('--backend', help='The backend to use. Can be one of local, spark, or batch')] = None,
           driver_cores: Ann[int, Opt('--driver-cores', help='Batch backend only. Number of cores to use for the driver process. May be 1, 2, 4, or 8.')] = 1,
           driver_memory: Ann[str, Opt('--driver-memory', help='Batch backend only. Memory tier to use for the driver process. May be standard or highmem.')] = 'standard',
           worker_cores: Ann[int, Opt('--worker-cores', help='Batch backend only. Number of cores to use for the worker processes. May be 1, 2, 4, or 8.')] = 1,
           worker_memory: Ann[str, Opt('--worker-memory', help='Batch backend only. Memory tier to use for the worker processes. May be standard or highmem.')] = 'standard',
           gcs_requester_pays_configuration: Ann[str, Opt('--gcs-requester-pays-configuration', help='If a string is provided, configure the Google Cloud Storage file system to bill usage to the project identified by that string. If a tuple is provided, configure the Google Cloud Storage file system to bill usage to the specified project for buckets specified in the list.')] = None,
           regions: Ann[List[str], Opt('--regions', help='List of regions to run jobs in when using the Batch backend. Use :data:`.ANY_REGION` to specify any region is allowed or use `None` to use the underlying default regions from the hailctl environment configuration. For example, use `hailctl config set batch/regions region1,region2` to set the default regions to use.')] = None,
           gcs_bucket_allow_list: Ann[str, Opt('--gcs-bucket-allow-list', help='A json list of buckets that Hail should be permitted to read from or write to, even if their default policy is to use "cold" storage.')] = None,
           copy_spark_log_on_error: Ann[bool, Opt('--copy-spark-log-on-error', help='Spark backend only. If `True`, copy the log from the spark driver node to `tmp_dir` on error.')] = False
           ):
    """Common Entry Point"""
    ctx.obj = HailInitKwargs(
        app_name=app_name,
        master=master,
        local=local,
        log=log,
        quiet=quiet,
        append=append,
        min_block_size=min_block_size,
        branching_factor=branching_factor,
        tmp_dir=tmp_dir,
        default_reference=default_reference,
        idempotent=idempotent,
        global_seed=global_seed,
        spark_conf=json.loads(spark_conf) if spark_conf else None,
        skip_logging_configuration=skip_logging_configuration,
        local_tmpdir=local_tmpdir,
        _optimizer_iterations=_optimizer_iterations,
        backend=backend,
        driver_cores=driver_cores,
        driver_memory=driver_memory,
        worker_cores=worker_cores,
        worker_memory=worker_memory,
        gcs_requester_pays_configuration=json.loads(gcs_requester_pays_configuration) if gcs_requester_pays_configuration else None,
        regions=regions,
        gcs_bucket_allow_list=json.loads(gcs_bucket_allow_list) if gcs_bucket_allow_list else None,
        copy_spark_log_on_error=copy_spark_log_on_error
    )


for cli in (
    qc_cli.qc_app,
):
    app.add_typer(cli)


def main():
    app()


if __name__ == '__main__':
    main()
