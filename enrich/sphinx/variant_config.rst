.. include:: global.rst

.. include:: seqlib_config.rst

**'wild type'** - *required*
	Information about the wild type sequence for this sequencing library.

	**'sequence'** - *required*
		The wild type DNA sequence. The sequence should be the same length as the sequences being compared to it (after any merging or read trimming).

	**'coding'** - *required*
		Set to **True** if the wild type sequence codes for protein. If this is **True**, amino acid changes will be computed for all variants. The wild type sequence must be in-frame.

	**'reference offset'**
		If this integer option is set, the value will be added to the DNA position (with respect to the wild type sequence) of each variant called. This is used to indicate the position of the mutagenized region within a longer sequence (mRNA, chromosome, etc.).

**'align variants'**
	Set to ``True`` to use the :py:class:`~seqlib.aligner.Aligner` to align variants that have too many mutations or an unexpected length. Calls indels as well as single nucleotide changes.

	.. note:: Alignment is typically disabled for performance reasons unless the user is interested in indel mutations.

