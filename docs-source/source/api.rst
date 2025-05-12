.. _sec-api:

====================
Python API Reference
====================

This is the API documentation for the BGE Toolkit, and provides detailed information
on the Python programming interface.

Use ``import bge_toolkit`` to access this functionality.


bge_toolkit.qc
~~~~~~~~~~~~~~

.. currentmodule:: bge_toolkit.qc

JointCallSet
============

A :class:`.JointCallSet` is an object that contains operations performed on the
combination of an exome Hail MatrixTable and an imputation Hail MatrixTable.

.. autoclass:: JointCallSet
    :members:
    :undoc-members:
    :inherited-members:
    :show-inheritance:


JoinType
========

An Enum that specifies what type of join to do. Inner joins are substantially cheaper
than outer joins!

.. autoclass:: JoinType
    :members:


Selectors
=========

.. data:: ALL_GROUP

   Use this global variable to select all variants and samples (no grouping).


.. data:: ALL_AGG

   Use this global variable to select all variants and samples for aggregation (no grouping).


Concordance
===========

A stand-alone Python function to compute concordance identical to the CLI

.. code-block::

    $ bge-toolit qc concordance

.. autofunction:: concordance


ConcordanceTable
================

A :class:`.ConcordanceTable` is an object that provides an interface for viewing concordance
results.

.. autoclass:: ConcordanceTable
    :members:
    :undoc-members:
    :inherited-members:
    :show-inheritance:


ConcordanceView
===============

A :class:`.ConcordanceView` is an object that provides an interface for viewing concordance
results for a specific grouping and set of aggregation variables.

.. autoclass:: ConcordanceView
    :members:
    :undoc-members:
    :inherited-members:
    :show-inheritance:


Statistic
=========

An Enum that specifies what type of concordance statistic to compute.

.. autoclass:: Statistic
    :members:


bge_toolkit.common
~~~~~~~~~~~~~~~~~~

.. currentmodule:: bge_toolkit.common

.. _binning-func:

Binning Functions
=================

.. autofunction:: mac_bin
.. autofunction:: maf_bin
.. autofunction:: gq_bin
.. autofunction:: dp_bin
.. autofunction:: max_gp_bin
.. autofunction:: qual_approx_bin
.. autofunction:: info_score_bin
