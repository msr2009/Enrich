.. include:: global.rst

Enrich2 Project To Do List
##########################

Debugging
=========

Output directory errors (from Doug)

Identify inconsistencies between old Enrich and new Enrich


Implementation - High Priority
==============================

Make a pass for adding additional status messages

Add more intelligent, complete output

* Write ``filter_stats`` files instead of writing to the log file

* Implement better :py:meth:`~datacontainer.DataContainer.write_data` based methods for the :py:class:`~selection.Selection` and :py:class:`~experiment.Experiment`

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

Implement wild-type correction ``[selection.py]``

Implement control-based correction ``[experiment.py]``

Implement filtering ``[experiment.py]``

Finish implementing "unlinking" ``[selection.py]``


Implementation - Low Priority
=============================

Extend logging to use multiple files (specifically to support sending filtered reads to their own output file)

Define a custom logging message formatting - "root" is unnecessary and confusing ``[enrich.py]``

Allow reading directly from ``.gz``/``.bz`` FASTQ_ files ``[fqread.py]``

Polish pass through ``import`` statements

Develop unit tests (use Python ``unittest`` library?)

Add multithreading to library read stage (one thread per SeqLib to speed up file I/O step)


Documentation
=============

Document logging behaviour (message types output at ``INFO``, ``DEBUG`` levels) and standard log message format::

	Capitalized message [self.name]

Write ``.rst`` file for the :py:class:`~datacontainer.DataContainer` class

Write a complete list of ``JSON`` configuration options for each level in a config file (each as its own ``.rst``?)

Documentation and ``.rst`` for ``trim_fastq.py``

Documentation and ``.rst`` for ``split_fastq.py``

Bare ``.html`` landing page for the documentation

More comprehensive description of filtering options
