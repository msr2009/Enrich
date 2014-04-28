.. include:: global.rst

.. include:: datacontainer_config.rst

**'timepoint'** - *required*
	Integer timepoint for this sequencing library. This can indicate time at which timepoints were taken (in hours, days, etc.) or the number of rounds of selection. The input library must be timepoint 0. Multiple sequencing libraries with the same timepoint will be combined by :py:class:`~selection.Selection`.

**'report filtered reads'**
	If this is ``True``, reads that are filtered out will be written to the log file, along with a message about which filters they failed.

	.. note:: Enabling this option can generate very large log files, and it is recommended that it should only be used for troubleshooting subsets of the data.

