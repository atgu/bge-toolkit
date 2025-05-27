Command-Line Interface
======================

bge-toolkit
-----------

.. custom-typer:: bge_toolkit.cli.__main__:app
    :prog: bge-toolkit

bge-toolkit qc
--------------

.. custom-typer:: bge_toolkit.qc.cli:qc_app
    :prog: bge-toolkit qc

bge-toolkit qc concordance
--------------------------

.. custom-typer:: bge_toolkit.qc.cli:concordance_app
    :prog: bge-toolkit qc concordance

============
Aggregations
============

``bge-toolkit qc concordance`` allows you to specify different grouping variables and aggregation variables when computing
global or variant concordance (not sample concordance).
Be aware the memory requirements / compute time / cost is proportional to the number of aggregations and grouping variables.
You can change the driver and work memory with the global flags ``--worker-memory`` and ``--driver-memory``.

The following aggregations are implemented:

Grouping
~~~~~~~~

- Exome MAF
- Exome MAC
- Imputation MAF
- Imputation MAC

Aggregation
~~~~~~~~~~~

- QualApprox
- INFO
- DP (computed as the sum of AD)
- GQ
- MAX_GP (computed as the maximum value of GP)


See :ref:`this section <binning-func>` for more details about the bins that are computed.

=====
Notes
=====

The following output files are written to ``--output-dir``:

1. Global Concordance

- global_concordance_table.ht (Contains the raw concordance data that can be loaded in Python with ``bge_toolkit.qc.ConcordanceTable.load()``)
- global-results/\*.tsv (Contains formatted tables for each combination of grouping and aggregation variable)
- global-plots/nonref-conc/\*.png (Contains figures of non-ref concordance for each grouping variable across aggregation variables)
- global-plots/f1-score/\*.png (Contains figures of f1-score for each grouping variable across aggregation variables)

2. Variant Concordance

- variant_overlaps.ht (Contains a Hail table with the list of overlapping variants)
- variant_conc.ht (Contains the raw concordance data that can be loaded in Python with ``bge_toolkit.qc.ConcordanceTable.load()``)

3. Sample Concordance

- sample_overlaps.ht (Contains a Hail table with the list of overlapping samples)
- sample_overlaps.tsv (Contains a TSV file with two columns for whether samples appear in the exome and imputation datasets)
- sample_conc.ht (Contains the raw concordance data that can be loaded in Python with ``bge_toolkit.qc.ConcordanceTable.load()``)
- sample-results/ALL_ALL.tsv (Contains the global concordance for each sample in the dataset)

========
Examples
========

1. Compute global concordance stratified by the Exome MAF with input files as MatrixTables for "chr1".

.. code-block::

    $ bge-toolkit \
         --driver-cores 2 \
         qc concordance \
         --exome "gs://my-bucket/exome.mt" \
         --imputation "gs://my-bucket/imputation.mt" \
         --output-dir "gs://my-bucket/concordance/" \
         --contig "chr1" \
         --EXOME-MAF \
         --global-conc


2. Compute global concordance stratified by Exome MAF and Imputation MAF with PLINK input files for a downsampled dataset.

.. code-block::

    $ bge-toolkit \
         --driver-cores 2 \
         qc concordance \
         --exome "gs://my-bucket/exome-plink" \
         --imputation "gs://my-bucket/imputation-plink" \
         --output-dir "gs://my-bucket/concordance/" \
         --downsample-variants 0.1 \
         --downsample-samples 0.1 \
         --EXOME-MAF \
         --IMPUTATION-MAF \
         --global-conc


3. Compute variant, sample, and global concordance statistics with VCF input files.

.. code-block::

    $ bge-toolkit \
         --driver-cores 2 \
         qc concordance \
         --exome "gs://my-bucket/exome.vcf.bgz" \
         --imputation "gs://my-bucket/imputation.vcf.bgz" \
         --output-dir "gs://my-bucket/concordance/" \
         --variant-conc \
         --sample-conc \
         --global-conc


bge-toolkit qc sample-qc
------------------------

.. custom-typer:: bge_toolkit.qc.cli:sample_qc_app
    :prog: bge-toolkit qc sample-qc

========
Examples
========

.. code-block::

    $ bge-toolkit qc sample-qc \
        --exome "gs://MY-BUCKET/data.mt" \
        --output-dir "gs://MY-BUCKET/test-sample-qc/080525-v2/" \
        --exome-regions "gs://MY-BUCKET/Twist_Alliance_Clinical_Research_Exome_Covered_Targets_hg38-34.9MB.bed" \
        --low-complexity-regions "gs://MY-BUCKET/LCRFromHengHg38.bed" \
        --dragen \
        --reported-sex-path "gs://MY-BUCKET/metadata.ht" \
        --reported-sex-col "reported_sex" \
        --chimera-rate-path "gs://MY-BUCKET/metadata.ht" \
        --chimera-rate-col "CHIMERA"

=====
Notes
=====

The following output files are written to ``--output-dir``:

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
