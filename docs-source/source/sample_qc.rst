.. _sec-sample_qc:

===================
Sample QC Reference
===================

This is the Python API documentation for the concordance component of the BGE Toolkit.

Use ``import bge_toolkit.qc`` to access this functionality.

.. currentmodule:: bge_toolkit.qc

sample_qc
---------

If you want to use the default implementation as described in Sealock et. al. "Tutorial: guidelines for quality filtering of whole-exome and whole-genome sequencing data for population-scale association analyses",
then you can use the function :func:`.sample_qc` with Hail MatrixTable and Table inputs.

An example of using the Python API is:

>>> from bge_toolkit.qc import sample_qc
>>> import hail as hl

>>> out_dir = "gs://my-bucket/my-sample-qc-output/"
>>> exome_regions_path = "gs://jigold-batch-tmp-ezxyx/Twist_Alliance_Clinical_Research_Exome_Covered_Targets_hg38-34.9MB.bed"
>>> low_complexity_regions_path = "gs://jigold-batch-tmp-ezxyx/LCRFromHengHg38.bed"

>>> exome_regions = hl.import_locus_intervals(exome_regions_path, reference_genome='GRCh38')
>>> low_complexity_regions = hl.import_locus_intervals(low_complexity_regions_path, reference_genome='GRCh38')
>>> exo = hl.read_matrix_table("gs://my-bucket/my-data.mt")

>>> sample_qc_stats = sample_qc(mt=exo, out_dir=out_dir, is_gatk=False, exome_regions=exome_regions, low_complexity_regions=low_complexity_regions)

>>> sample_qc_stats.write(f'{out_dir}sample_qc.ht')
>>> sample_qc_stats.pcs.export(f'{out_dir}sample_pcs.tsv')
>>> sample_qc_stats.plot_pcs(1, 2)
>>> sample_qc_stats.plot_gq_fraction_density()
>>> sample_qc_stats.plot_qc_metric_boxplots()['r_ti_tv']

The following output files are written to ``out_dir``:

- sample_qc_stats.ht (Contains the raw sample qc data that can be loaded in Python with :meth:`.SampleQCResult.load`)
- pcs/pc1_pc2.png (Contains a plot of PC1 versus PC2 colored by ancestry population label)
- pcs/pc1_pc3.png (Contains a plot of PC1 versus PC3 colored by ancestry population label)
- pcs/pc2_pc3.png (Contains a plot of PC2 versus PC3 colored by ancestry population label)
- qc/\*_boxplot.png (Contains boxplots of different QC metrics stratified by ancestry population label with outliers flagged)
- qc/\*_density.png (Contains density plots of different QC metrics stratified by ancestry population label)
- passing_sample_ids.tsv (A TSV file containing a list of sample IDs for samples that passed all QC metrics)
- qc/\*_pass_boxplot.png (Contains boxplots of different QC metrics stratified by ancestry population label with outliers flagged for passing samples only)
- qc/\*_pass_density.png (Contains density plots of different QC metrics stratified by ancestry population label for passing samples only)

The structure of ``sample_qc_stats.ht`` is as follows:

.. code-block::

    contamination:
        - charr: The CHARR statistic.
        - is_passing: ``charr`` is less than the ``charr_thresh``

    chimera_reads:
        - chimera_rate: The rate of chimera reads.
        - is_passing: The rate of chimera reads is below ``threshold``.

    sample_qc_metrics:
        - is_passing: Passes every statistic.
        - r_ti_tv: Transition / Transversion ratio.
        - n_singleton: Number of singletons.
        - n_insertion: Number of insertions.
        - n_deletion: Number of deletions.
        - n_transition: Number of transitions.
        - n_transversion: Number of transversions.
        - r_het_hom_var: Ratio of heterozygotes to number of homozygote variants.
        - r_insertion_deletion: Ratio of insertions to deletions.
        - fail_r_ti_tv: An outlier in ``r_ti_tv``.
        - fail_n_singleton: An outlier in ``n_singleton``.
        - fail_n_insertion: An outlier in ``n_insertion``.
        - fail_n_deletion: An outlier in ``n_deletion``.
        - fail_n_transition: An outlier in ``n_transition``.
        - fail_n_transversion: An outlier in ``n_transversion``.
        - fail_r_het_hom_var: An outlier in ``r_het_hom_var``.
        - fail_r_insertion_deletion: An outlier in ``r_insertion_deletion``.

    gq_fraction:
        - gq_fraction: The percentage of genotypes with GQ >= ``gq_thresh``.
        - is_passing: ``gq_fraction`` is greater than the ``fraction_thresh``

    sex_info:
        - is_female: The imputed sex. ``True`` is for females, ``False`` is for males.
        - is_passing: Either ``reported sex == imputed sex`` or ``True``.
        - sex_check: ``reported sex == imputed sex`` or ``Null``.
        - f_stat: Inbreeding coefficient on the non-PAR X chromosome.
        - n_called: Number of genotypes considered.
        - expected_homs: Expected number of homozygotes.
        - observed_homs: Observed number of heterozygotes.

    ancestry:
        - ancestry_pop: The ancestry population label.

    pcs:
        - scores: An array with the top ``k`` principal components.

    s: sample ID


.. autofunction:: sample_qc


High Level Functions
--------------------

The :func:`.sample_qc` function is implemented in terms of the following three lower-level functions:

.. autofunction:: select_high_quality_common_sites
.. autofunction:: infer_ancestry
.. autofunction:: calculate_sample_qc_stats


Lower Level Functions
---------------------

:func:`.select_high_quality_common_sites`, :func:`.infer_ancestry`, and :func:`.calculate_sample_qc_stats`
are implemented in terms of these lower level functions which you can use to build your own QC pipeline:

.. autofunction:: filter_to_high_quality_exome_intervals
.. autofunction:: filter_contamination_rate
.. autofunction:: filter_chimeric_reads
.. autofunction:: filter_genotypes
.. autofunction:: filter_to_high_quality_common_sites
.. autofunction:: compute_pcs
.. autofunction:: estimate_call_rate_through_gq_stats
.. autofunction:: impute_and_check_sex
.. autofunction:: identify_outliers_in_qc_stats


SampleQCResult
--------------

A :class:`.SampleQCResult` is an object that provides a wrapper for interacting
with the sample QC results by making plots and exploring the metrics.

.. autoclass:: SampleQCResult
    :members:
    :undoc-members:
    :inherited-members:
    :show-inheritance:
