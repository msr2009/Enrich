.. include:: global.rst

Enrich2 Project To Do List
##########################

Debugging
=========

Write unit tests for all files

* :py:class:`~fqread.FQRead`
* trim_fastq.py
* split_fastq.py
* config_check.py
* enrich_error.py
* :py:class:`~datacontainer.DataContainer`
* :py:class:`~seqlib.seqlib.SeqLib`
* :py:class:`~seqlib.variant.VariantSeqLib`
* :py:class:`~seqlib.barcode.BarcodeSeqLib`
* :py:class:`~seqlib.basic.BasicSeqLib`
* :py:class:`~seqlib.overlap.OverlapSeqLib`
* :py:class:`~seqlib.barcodevariant.BarcodeVariantSeqLib`
* :py:class:`~selection.Selection`
* :py:class:`~experiment.Experiment`

Implementation - High Priority
==============================

Make a pass for adding additional status messages

Add basic plotting functions

* :py:class:`~seqlib.seqlib.SeqLib`
	* Diversity heatmap
	* Frequency histogram

* :py:class:`~selection.Selection`
	* Sequence-function map
	* Score histogram

* :py:class:`~experiment.Experiment`
	* Sequence-function map

Implement replicate comparison ``[experiment.py]``

Implement control-based correction ``[experiment.py]``

Implement filtering ``[experiment.py]``

Finish implementing "unlinking" ``[selection.py]``


Implementation - Low Priority
=============================

Extend logging to use multiple files (specifically to support sending filtered reads to their own output file)

Define a custom logging message formatting - "root" is unnecessary and confusing ``[enrich.py]``

Polish pass through ``import`` statements

Add multithreading to library read stage (one thread per SeqLib to speed up file I/O step)

Move to Python3-style (.format) string formatting, as is recommended for new code


Documentation
=============

Document logging behaviour (message types output at ``INFO``, ``DEBUG`` levels) and standard log message format::

	Capitalized message [self.name]

Bare ``.html`` landing page for the documentation

More comprehensive description of filtering options

Better description of the output file structure (especially *subdirectory* for :py:meth:`~datacontainer.DataContainer.write_data`)
