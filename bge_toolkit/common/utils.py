import logging
import os
from typing import List, Optional, Tuple, Union

import hail as hl
import hailtop.fs as hfs
from hail.backend.service_backend import ServiceBackend
from plotnine import ggplot
from tempfile import NamedTemporaryFile


def setup_logger(
    log_name: str = __name__,  # Use __name__ for a unique, standard logger name
    log_path: Optional[str] = None,
    level: int = logging.DEBUG
) -> logging.Logger:
    """
    Sets up a logger with a stream handler and an optional file handler.
    This function is idempotent and safe for use in a library.
    """
    log = logging.getLogger(log_name)

    # Only configure if the logger hasn't been configured yet
    if not log.hasHandlers():
        log.setLevel(level)
        log.propagate = False

        # Create and add the stream handler for console output
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(level)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        stream_handler.setFormatter(formatter)
        log.addHandler(stream_handler)

        # Optionally create and add the file handler
        if log_path:
            file_handler = logging.FileHandler(log_path, mode='w')
            file_handler.setLevel(level)
            file_handler.setFormatter(formatter)
            log.addHandler(file_handler)

    return log


def save_ggplot(p: ggplot, path: str, **kwargs):
    if path is not None:
        with NamedTemporaryFile(delete=False, suffix='.png') as temp_file:
            temp_file_name = temp_file.name

        p.draw()
        p.save(temp_file_name, transparent=False, **kwargs)

        hfs.copy(temp_file_name, path)
        os.remove(temp_file_name)


def setup_logging_and_qob(log_name: str,
                          log_path: Optional[str] = None,
                          hail_init_kwargs: Optional[dict] = None,
                          log: Optional[logging.Logger] = None) -> logging.Logger:
    log = setup_logger(log, log_path, log_name)
    if log is None:
        log = logging.getLogger(log_name)
        log.setLevel(logging.INFO)
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


def apply_filters(*,
                  dataset: Union[hl.MatrixTable, hl.vds.VariantDataset],
                  description: Optional[str] = None,
                  log: Optional[logging.Logger] = None,
                  sample_list: Optional[str] = None,
                  variant_list: Optional[str] = None,
                  contig: Optional[List[str]] = None,
                  n_samples: Optional[int] = None,
                  n_variants: Optional[int] = None,
                  downsample_samples: Optional[float] = None,
                  downsample_variants: Optional[float] = None) -> Tuple[hl.MatrixTable, bool]:
    """Apply filters to MatrixTable.

    Args:
        dataset (MatrixTable, VariantDataset): MatrixTable or VariantDataset to filter
        description (Optional[str]): Description of MatrixTable such as the original path.
        log (Optional[logging.Logger]): A logging object to log to.
        sample_list (Optional[str]): Filter to samples listed in the file.
        variant_list (Optional[str]): Filter to variants listed in the file.
        contig (Optional[List[str]]): Filter to variants in the contig.
        n_samples (Optional[int]): Filter to the first N overlapping samples.
        n_variants (Optional[int]): Filter to the first N variants.
        downsample_samples (Optional[float]): Downsample to X fraction of samples.
        downsample_variants (Optional[float]): Downsample to X fraction of variants.

    Returns:
        A filtered MatrixTable or VariantDataset along with a boolean specifying whether the dataset has been downsampled.
    """
    downsampled = False

    description = description or '<Unknown Source>'

    if log is None:
        log = setup_logger(None, None, 'filter_dataset')

    if sample_list is not None:
        samples = hl.import_table(sample_list)
        samples = samples.key_by('s')
        if isinstance(dataset, hl.MatrixTable):
            dataset = dataset.semi_join_cols(samples)
        else:
            assert isinstance(dataset, hl.vds.VariantDataset)
            dataset = hl.vds.filter_samples(dataset, samples, keep=True)
        log.info(f'Subset to samples in {sample_list} for {description}.')

    if variant_list is not None:
        downsampled = True
        variants = hl.import_table(variant_list, no_header=True)
        variants = variants.key_by('v')
        variants = variants.annotate(variant=hl.parse_variant(variants.v))
        variants = variants.annotate(locus=variants.variant.locus, alleles=variants.variant.alleles)
        variants = variants.key_by('locus', 'alleles')

        if isinstance(dataset, hl.MatrixTable):
            dataset = dataset.semi_join_rows(variants)
        else:
            assert isinstance(dataset, hl.vds.VariantDataset)
            variant_data = dataset.variant_data
            variant_data = variant_data.semi_join_rows(variants)
            dataset = hl.vds.VariantDataset(dataset.reference_data, variant_data)

        log.info(f'Subset to variants in {variant_list} for {description}.')

    if contig is not None:
        downsampled = True
        if isinstance(dataset, hl.MatrixTable):
            contig_expr = dataset.locus.contig == contig[0]
            for c in contig[1:]:
                contig_expr |= dataset.locus.contig == c
            dataset = dataset.filter_rows(contig_expr)
        else:
            assert isinstance(dataset, hl.vds.VariantDataset)
            dataset = hl.vds.filter_chromosomes(dataset, keep=contig)
        log.info(f'Subsetting to variants on {contig} for {description}.')

    if downsample_variants:
        downsampled = True
        if isinstance(dataset, hl.MatrixTable):
            dataset = dataset.sample_rows(downsample_variants, seed=10243523)
        else:
            assert isinstance(dataset, hl.vds.VariantDataset)
            variant_data = dataset.variant_data.sample_rows(downsample_variants, seed=10243523)
            dataset = hl.vds.VariantDataset(dataset.reference_data, variant_data)
        log.info(f'Downsampling variants by {downsample_variants} for {description}.')

    if downsample_samples:
        downsampled = True
        if isinstance(dataset, hl.MatrixTable):
            dataset = dataset.sample_cols(downsample_samples, seed=10243523)
        else:
            assert isinstance(dataset, hl.vds.VariantDataset)
            samples_to_downsample = dataset.variant_data.sample_cols(downsample_samples, seed=10243523).cols()
            dataset = hl.vds.filter_samples(dataset, samples_to_downsample)
        log.info(f'Downsampling samples by {downsample_samples} for {description}.')

    if n_variants is not None:
        downsampled = True
        if isinstance(dataset, hl.MatrixTable):
            dataset = dataset.head(n_rows=n_variants)
        else:
            assert isinstance(dataset, hl.vds.VariantDataset)
            variant_data = dataset.variant_data.head(n_rows=n_variants)
            dataset = hl.vds.VariantDataset(dataset.reference_data, variant_data)
        log.info(f'Selecting the first {n_variants} variants for {description}.')

    if n_samples is not None:
        downsampled = True
        if isinstance(dataset, hl.MatrixTable):
            dataset = dataset.head(n_rows=None, n_cols=n_samples)
        else:
            assert isinstance(dataset, hl.vds.VariantDataset)
            samples_to_keep = dataset.variant_data.head(n_rows=None, n_cols=n_samples).cols()
            dataset = hl.vds.filter_samples(dataset, samples_to_keep, keep=True)
        log.info(f'Selecting the first {n_samples} samples for {description}.')

    return (dataset, downsampled)
