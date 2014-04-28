.. include:: global.rst

.. include:: variant_config.rst

**'fastq'** - *required*
	Information about the FASTQ_ file.

	**'forward'** and **'reverse'** - *required*
		Both FASTQ_ files must be specified. All reads in the 'reverse' file will be reverse-complemented during the merging process.

**'overlap'** - *required*
	Information about how the forward and reverse reads should be combined.

	**'forward start'** - *required*
		Position in the forward read (first base is position 1) where the overlapping region begins.

	**'reverse start'** - *required*
		Position in the reverse read (first base is position 1) where the overlapping region begins. This is the position in the read before reverse-complementing.

	**'length'** - *required*
		Number of bases in the overlaping region.

**'overlap only'** - *required*
	If this option is ``True``, the portions of the read that are outside the overlapping region will be trimmed before calling variants. Otherwise, the portion of the forward read before the overlapping region and the portion of the reverse read after the overlapping region are included in the variant sequence.

**'max mismatches'** - *required*
	The maximum number of mismatches tolerated in the overlapping region.

**'filters'** - *required*
	Filtering options for reads and variants.

	**'min quality'**
		Minimum quality value for all bases in the merged read.

	**'avg quality'**
		Minimum average quality value for the merged read.

	**'chastity'**
		If ``True``, require that both forward and reverse reads have the chastity bit set in their FASTQ_ headers.

	**'max mutations'**
		Maximum number of mutations allowed for the variant.

	**'remove unresolvable'**
		Remove merged reads with unresolvable mismatches (different nucleotides with the same quality score at the same position).


