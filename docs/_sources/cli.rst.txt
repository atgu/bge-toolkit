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
