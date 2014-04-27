.. include:: global.rst

:py:mod:`~seqlib.aligner` --- Needleman-Wunsch alignment for variants
=====================================================================

.. py:module:: seqlib.aligner
    :synopsis: Needleman-Wunsch alignment for variants.

The :py:mod:`~seqlib.aligner` module contains a `Needleman-Wunsch 
<http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm>`_ aligner 
used by :py:class:`~seqlib.variant.VariantSeqLib` objects to align variants 
to the wild type sequence.

.. note:: This module uses Python's :py:class:`Exception` instead of :py:class:`~enrich_error.EnrichError` for portability.

:py:class:`~seqlib.aligner.Aligner` class
-----------------------------------------
.. autoclass:: Aligner(similarity=_simple_similarity)
    :members:
