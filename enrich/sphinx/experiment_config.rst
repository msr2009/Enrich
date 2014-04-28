.. include:: global.rst

.. include:: datacontainer_config.rst

**'conditions'** - *required*
	List of experimental conditions, each having the following format.

	**'label'** - *required*
		Label for the condition. Must be unique and alphanumeric.

	**'selection'** - *required*
		List of :py:class:`~selection.Selection` config objects for selection performed under this condition.

	**'control'**
		If this is ``True``, this condition will be treated as the control for normalization purposes. Only one control condition may be specified.

**'normalize wt'**
	If this is ``True``, normalize all scores or ratios using the wild type score or ratio as neutral.
