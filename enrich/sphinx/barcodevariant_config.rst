.. include:: global.rst

.. include:: variant_config.rst

**'fastq'** - *required*
	Information about the FASTQ_ file.

	**'forward'** or **'reverse'** - *required*
		Only one FASTQ_ file may be specified. If the file is 'reverse', all reads will be reverse-complemented before barcodes are counted.

	**'start'**
		Starting position of the barcode in the read (first base is position 1). Used for optional read trimming.

	**'length'**
		Number of bases in the barcode. Used for optional read trimming.

**'barcodes'** - *required*
	This config option must be present for the sequences to be treated as barcodes, even if it has no elements in it.

	**'map file'**
		Path to the :py:class:`~seqlib.barcodevariant.BarcodeMap` file mapping barcodes to variant sequences. If this option is not set, the :py:class:`~seqlib.barcodevariant.BarcodeMap` from the :py:class:`~selection.Selection` will be used.

	**'min count'**
		Minimum count for a barcode to be included in the analysis. Barcodes with counts below this threshold will be output as low abundance barcodes, then discarded.

**'filters'** - *required*
	Filtering options for reads and variants.

	**'min quality'**
		Minimum quality value for all bases in the read.

	**'avg quality'**
		Minimum average quality value for the read.

	**'chastity'**
		If ``True``, require that the read has the chastity bit set in the FASTQ_ header.

	**'max mutations'**
		Maximum number of mutations allowed for the variant.

