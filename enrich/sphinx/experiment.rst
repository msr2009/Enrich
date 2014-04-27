.. include:: global.rst

:py:mod:`~experiment` --- Class for comparing multiple selections
=================================================================

.. py:module:: experiment
	:synopsis: Class for comparing multiple selections.

The :py:mod:`~experiment` module contains the class definition for the :py:class:`~experiment.Experiment` class, used for comparing multiple :py:class:`~selection.Selection` objects. Typically, one :py:class:`~experiment.Experiment` will be created to perform all analysis for a single run of the `Enrich2 <index.html>`_ pipeline.

Configuration options
---------------------

.. include:: experiment_config.rst

:py:class:`~experiment.Experiment` class
----------------------------------------
.. autoclass:: experiment.Experiment
	:members:

:py:mod:`~experiment` apply functions
-------------------------------------
.. autofunction:: condition_cv_apply_fn
