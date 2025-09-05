import abc
import os
from typing import Optional
import logging

import hail as hl
import hailtop.fs as hfs
from hail.ggplot.utils import n_partitions
from hailtop.fs.fs_utils import GCSRequesterPaysConfiguration

log = logging.getLogger('concordance')


class OutputData(abc.ABC):
    typ: str

    def export(self, mt: hl.MatrixTable):
        pass


class VCFOutputData(OutputData):
    typ = 'vcf'

    def __init__(self, output_path: str, **kwargs):
        self.output_path = output_path
        self.kwargs = kwargs

    def export(self, mt: hl.MatrixTable):
        output_path = f'{self.output_path}.vcf.bgz'
        hl.export_vcf(mt, output_path, tabix=True, **self.kwargs)


class PlinkOutputData(OutputData):
    typ = 'plink'

    def __init__(self, output_path: str, **kwargs):
        self.output_path = output_path
        self.kwargs = kwargs

    def export(self, mt: hl.MatrixTable,):
        output_path = f'{self.output_path}.vcf.bgz'
        hl.export_plink(mt, output_path, ind_id=mt.col_key[0])


class MatrixTableOutputData(OutputData):
    typ = 'hail_mt'

    def __init__(self, output_path: str, overwrite: bool = False):
        self.output_path = output_path
        self.overwrite = overwrite

    def export(self, mt: hl.MatrixTable):
        mt.write(self.output_path, overwrite=self.overwrite)


class InputData(abc.ABC):
    typ: str
    _dataset: Optional[Union[hl.MatrixTable, hl.vds.VariantDataset]]

    def load(self) -> hl.MatrixTable:
        pass

    def count(self) -> Tuple[int, int]:
        pass


class PlinkInputData(InputData):
    typ = 'plink'

    def __init__(self, bfile_root, **kwargs):
        self.bfile_root = bfile_root
        self.kwargs = kwargs
        self._dataset = None

    def load(self) -> hl.MatrixTable:
        if self._dataset is None:
            mt = hl.import_plink(bed=f'{self.bfile_root}.bed',
                                 bim=f'{self.bfile_root}.bim',
                                 fam=f'{self.bfile_root}.fam',
                                 **self.kwargs)
            self._dataset = mt
        return self._dataset

    def count(self) -> Tuple[int, int]:
        data = self.load()
        return data.count()


class VCFInputData(InputData):
    typ = 'vcf'

    def __init__(self, path: str, **kwargs):
        self.path = path
        self.kwargs = kwargs
        self._dataset = None

    def load(self) -> hl.MatrixTable:
        if self._dataset is None:
            if self.path[-4:] == '.bgz' or self.path[-3:] == '.gz':
                data = hl.import_vcf(self.path,
                                     force_bgz=True,
                                     call_fields=['LGT'],
                                     **self.kwargs)
            else:
                data = hl.import_vcf(self.path, call_fields=['LGT'], **self.kwargs)

            if 'GT' not in data.entry and 'LGT' in data.entry and 'LA' in data.entry:
                data = data.annotate_entries(
                    GT=hl.vds.lgt_to_gt(data.LGT, data.LA)
                )

            self._dataset = data

        return self._dataset

    def count(self) -> Tuple[int, int]:
        data = self.load()
        return data.count()


class MatrixTableInputData(InputData):
    typ = 'hail_mt'

    def __init__(self, path: str, n_partitions: Optional[int]):
        self.path = path
        self.n_partitions = n_partitions
        self._dataset = None

    def load(self) -> hl.MatrixTable:
        if self._dataset is None:
            data = hl.read_matrix_table(self.path, _n_partitions=self.n_partitions)

            if 'GT' not in data.entry and 'LGT' in data.entry and 'LA' in data.entry:
                data = data.annotate_entries(
                    GT=hl.vds.lgt_to_gt(data.LGT, data.LA)
                )

            self._dataset = data

        return self._dataset

    def count(self) -> Tuple[int, int]:
        data = self.load()
        return data.count()


class VariantDatasetInputData(InputData):
    typ = 'hail_vds'

    def __init__(self, path: str, n_partitions: Optional[int]):
        self.path = path
        self.n_partitions = n_partitions
        self._dataset = None

    def load(self) -> hl.VariantDataset:
        if self._dataset is None:
            data = hl.vds.read_vds(self.path, n_partitions=self.n_partitions)

            variant_data = data.variant_data

            if 'GT' not in variant_data.entry and 'LGT' in variant_data.entry and 'LA' in variant_data.entry:
                variant_data = variant_data.annotate_entries(
                    GT=hl.vds.lgt_to_gt(variant_data.LGT, variant_data.LA)
                )

            data = hl.vds.VariantDataset(data.reference_data, variant_data)

            self._dataset = data

        return self._dataset

    def count(self) -> Tuple[int, int]:
        data = self.load()
        return data.variant_data.count()


class SplitRowsMatrixTableInputData(InputData):
    typ = 'split_rows_hail_mt'

    def __init__(self, path: str, requester_pays_config: Optional[GCSRequesterPaysConfiguration] = None):
        self.path = path
        self.requester_pays_config = requester_pays_config
        self._dataset = None

    def load(self) -> hl.MatrixTable:
        if self._dataset is None:
            file_entries = hfs.ls(self.path, requester_pays_config=self.requester_pays_config)
            paths = [entry.path for entry in file_entries if entry.path.endswith('.mt')]

            mt = hl.read_matrix_table(paths[0])
            datasets = [hl.read_matrix_table(p) for p in paths[1:]]

            data = mt.union_rows(*datasets)

            if 'GT' not in data.entry and 'LGT' in data.entry and 'LA' in data.entry:
                data = data.annotate_entries(
                    GT=hl.vds.lgt_to_gt(data.LGT, data.LA)
                )

            self._dataset = data

        return self._dataset

    def count(self) -> Tuple[int, int]:
        data = self.load()
        return data.count()


def load_input_data(path: str, **kwargs) -> InputData:
    if path[-8:] == '.vcf.bgz' or path[-7:] == '.vcf.gz' or path[-4:] == '.vcf':
        log.info(f'Found VCF input file: {path}')
        return VCFInputData(path, **kwargs)

    if path[-3:] == '.mt':
        log.info(f'Found MatrixTable input file: {path}')
        if '*' in path:
            return SplitRowsMatrixTableInputData(os.path.dirname(path)).load()
        return MatrixTableInputData(path, n_partitions=kwargs.get('n_partitions'))

    if path[-3:] == '.vds':
        log.info(f'Found VariantDataset input file: {path}')
        return VariantDatasetInputData(path)

    assert hfs.exists(f'{path}.bed') and hfs.exists(f'{path}.bim') and hfs.exists(f'{path}.fam')
    log.info(f'Found PLINK input files: {path}')
    return PlinkInputData(path)


def load_data(*, path: str, descriptor: str, log: logging.Logger, n_partitions: Optional[int] = None) -> hl.MatrixTable:
    log.info(f'loading dataset {descriptor}')
    input = load_input_data(path, n_partitions=n_partitions)
    n_variants, n_samples = input.count()
    log.info(f'found {n_variants} variants and {n_samples} samples in {descriptor}.')
    return input.load()
