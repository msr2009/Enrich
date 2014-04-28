.. include:: global.rst

**'name'**  - *required*
	The name of this entity. This used for output directory names, and is included in every error message and log message. Names should be unique within a single analysis, but this is not enforced.

**'output directory'**
	Sets the base output directory for the analysis. Output files and directories will be created in subdirectories at this location. It will be created if it doesn't exist.

	If desired, the output directory can be set for each entity separately. If no 'output directory' is included in the config, the entity will use the parent's setting.

	This setting is required for the top level entity (usually an :py:class:`~experiment.Experiment`) in each analysis.

