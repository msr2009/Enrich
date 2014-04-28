.. include:: global.rst

.. include:: variant_config.rst

**'fastq'** - *required*
	Information about the FASTQ_ file.

	**'forward'** or **'reverse'** - *required*
		Only one FASTQ_ file may be specified. If the file is 'reverse', all reads will be reverse-complemented before variants are called.

**'filters'**  *required*
	Filtering options for reads and variants.

	**'min quality'**
		Minimum quality value for all bases in the read.

	**'avg quality'**
		Minimum average quality value for the read.

	**'chastity'**
		If ``True``, require that the read has the chastity bit set in the FASTQ_ header.

	**'max mutations'**
		Maximum number of mutations allowed for the variant.

