.. include:: global.rst

:py:mod:`~selection` --- Class for analyzing a selection timeseries
===================================================================

.. py:module:: selection
	:synopsis: Class for analyzing a selection timeseries.

The :py:mod:`~selection` module contains the class definition for the :py:class:`~selection.Selection` class, used for analyzing multiple timepoints in a selection experiment. Each timepoint consists of one or more `Seqlib family objects <seqlib_overview.html>`_.

Configuration options
---------------------

.. include:: selection_config.rst

:py:class:`~selection.Selection` class
--------------------------------------
.. autoclass:: selection.Selection
	:members:

:py:mod:`~selection` apply functions
-------------------------------------
.. autofunction:: nonsense_ns_carryover_apply_fn

.. autofunction:: enrichment_apply_fn

Filtering apply functions
*************************
.. autofunction:: min_count_filter

.. autofunction:: min_input_count_filter

.. autofunction:: min_rsq_filter

:py:class:`~seqlib.barcodevariant.BarcodeVariant`-specific apply functions
**************************************************************************
.. autofunction:: barcode_variation_apply_fn

.. autofunction:: barcode_count_apply_fn

.. autofunction:: barcode_varation_filter
