.. include:: global.rst

.. include:: datacontainer_config.rst

**'libraries'** - *required*
	List of :py:class:`~seqlib.seqlib.SeqLib` config objects for sequencing libraries that are part of this :py:class:`~selection.Selection`.

**'barcodes'**
	Barcode options.

	**'map file'**
		Path to the :py:class:`~seqlib.barcodevariant.BarcodeMap` file mapping barcodes to variant sequences. This :py:class:`~seqlib.barcodevariant.BarcodeMap` will be used for any :py:class:`~seqlib.barcodevariant.BarcodeVariantSeqLib` libraries that do not have their own specified.

**'filters'** - *required*
	Filtering options for barcodes and variants.

	**'min count'**
		Minimum count for this barcode or variant at any timepoint.

	**'min input count'**
		Minimum count for this barcode or variant in the input (timepoint 0).

	**'min rsquared'**
		Minimum r-squared value for scored variants and barcodes.

	**'max barcode variation'**
		Maximum coefficient of variation for barcode scores mapped to the same variant. This is only applied when only when all sequencing data has both barcode and variant data (*i.e.* are :py:class:`~seqlib.barcodevariant.BarcodeVariantSeqLib` objects) with the same :py:class:`~seqlib.barcodevariant.BarcodeMap`.

**'carryover correction'**
	Nonspecific carryover correction options.

	**'method'**
		Nonspecific carryover correction method. Currently, only "nonsense" is supported (described in `Araya and Fowler`_).

	**'position'**
		Additional parameterd for "nonsense" method. This is the last amino acid position (first position is 1) for which a change to stop marks a variant as part of nonspecific carryover.

**'normalize wt'**
	If this is ``True``, normalize all scores or ratios using the wild type score or ratio as neutral.
