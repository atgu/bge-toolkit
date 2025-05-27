.. _sec-api:

====================
Python API Reference
====================

This is the API documentation for the BGE Toolkit, and provides detailed information
on the Python programming interface.

Use ``import bge_toolkit`` to access this functionality.


bge_toolkit.qc
--------------

.. currentmodule:: bge_toolkit.qc

The concordance functionality in the BGE Toolkit is meant to compute
concordance between an exome sequencing dataset and the corresponding imputed
dataset across arbitrary metrics to identify appropriate variant QC thresholds.


.. toctree::
    :maxdepth: 4

    concordance


The Sample QC functionality in the BGE Toolkit is meant to run standard
analyses for computing sample QC statistics. The pipeline is broken into
three components:

1. Find high quality variants.
2. Use a random forest model from gnomAD to identify ancestry population labels per sample.
3. Compute summary statistics including contamination rates, principal components, sex check, estimating call rate via GQ, and running `hl.sample_qc()` and identifying outliers by ancestry population.

Note: Relatedness checks have not been implemented yet due to a memory leak in Hail QoB.


.. toctree::
    :maxdepth: 4

    sample_qc
