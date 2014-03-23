from __future__ import print_function
import re
import logging
from variant import VariantSeqLib
from barcode import BarcodeSeqLib
from seqlib import SeqLib
from enrich_error import EnrichError
from fqread import read_fastq, check_fastq
import pandas as pd


class BarcodeMap(dict):
    """
    Dictionary-derived class for storing the relationship between barcodes 
    and variants. Requires the path to a *mapfile*, containing lines in the 
    format ``'barcode<tab>variant'`` for each barcode expected in the library. 
    Also creates a second dictionary, ``BarcodeMap.variants``, storing a list 
    of barcodes assigned to a given variant.

    Barcodes must only contain the characters ``ACGT`` and variants must only 
    contain the characters ``ACGTN`` (lowercase characters are also accepted). 

    The dictionaries are created when the object is initialized.
    """
    def __init__(self, mapfile):
        self.name = "mapfile_%s" % mapfile
        try:
            handle = open(mapfile, "U")
        except IOError:
            raise EnrichError("Could not open barcode map file '%s'" \
                    % mapfile, self.name)

        self.filename = mapfile
        for line in handle:
            # skip comments and whitespace-only lines
            if len(line.strip()) == 0 or line[0] == '#':
                continue

            try:
                barcode, variant = line.strip().split()
            except ValueError:
                raise EnrichError("Unexpected barcode-variant line format", 
                                  self.name)

            if not re.match("^[ACGTacgt]+$", barcode):
                raise EnrichError("Barcode DNA sequence contains unexpected "
                                  "characters", self.name)
            if not re.match("^[ACGTNacgtn]+$", variant):
                raise EnrichError("Variant DNA sequence contains unexpected "
                                  "characters", self.name)

            barcode = barcode.upper()
            variant = variant.upper()
            if barcode in self:
                if self[barcode] != variant:
                    raise EnrichError("Barcode '%s' assigned to multiple "
                                      "unique variants" % barcode, self.name)
            else:
                self[barcode] = variant
        handle.close()

        # build the variants dictionary
        self.variants = dict()
        for bc in self.keys():
            if self[bc] not in self.variants:
                self.variants[self[bc]] = list()
            self.variants[self[bc]].append(bc)

        logging.info("Assigned %d barcodes to %d variants [%s]" % \
                     (len(self.keys()), len(self.variants.keys()), self.name))


    def write_variants(self, fname):
        """
        Write a list of barcodes for each variant to the file *fname*.
        """
        try:
            handle = open(fname, "w")
        except IOError:
            raise EnrichError("Could not open variant barcode map file '%s' "
                              "for writing" % fname, self.name)
        for variant, barcodes in \
                sorted(self.variants.items(), key=lambda x:x[1]):
            print(variant, ", ".join(barcodes), sep="\t", file=handle)
        handle.close()

        logging.info('Wrote BarcodeMap variants file "%s" [%s]' % (fname, self.name))



class BarcodeVariantSeqLib(VariantSeqLib, BarcodeSeqLib):
    """
    Class for counting variant data from barcoded sequencing libraries. 
    Creating a :py:class:`BarcodeVariantSeqLib` requires a valid *config* 
    object with an ``'barcodes'`` entry and information about the wild type 
    sequence.

    Example config file for a :py:class:`~seqlib.barcodevariant.BarcodeVariantSeqLib`:

    .. literalinclude:: config_examples/barcodevariant.json

    :download:`Download this JSON file <config_examples/barcodevariant.json>`

    The ``barcode_map`` keyword argument can be used to pass an existing 
    :py:class:`~seqlib.barcodevariant.BarcodeMap`, but only if the 
    ``"map file"`` entry is absent from *config*.
    """
    def __init__(self, config, barcode_map=None):
        VariantSeqLib.__init__(self, config)
        BarcodeSeqLib.__init__(self, config, barcodevariant=True)
        try:
            if 'map file' in config['barcodes']:
                self.barcode_map = BarcodeMap(config['barcodes']['map file'])
            else:
                self.barcode_map = None

            self.set_filters(config['filters'], {'min quality' : 0,
                                      'avg quality' : 0,
                                      'chastity' : False,
                                      'max mutations' : len(self.wt_dna)})
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, 
                              self.name)

        if self.barcode_map is None: # not in local config
            if barcode_map is None:  # not provided on object creation
                raise EnrichError("Barcode map not specified", self.name)
            else:
                self.barcode_map = barcode_map

        self.df_dict['barcodes_unmapped'] = None
        self.filter_unmapped = True


    def calculate(self):
        """
        Counts the barcodes using :py:meth:`BarcodeSeqLib.count` and combines them into 
        variant counts using the :py:class:`BarcodeMap`.
        """
        BarcodeSeqLib.calculate(self) # count the barcodes
        self.df_dict['variants'] = dict()

        logging.info("Converting barcodes to variants [%s]" % self.name)
        if self.filter_unmapped:
            map_mask = self.df_dict['barcodes'].index.isin(self.barcode_map)
            self.df_dict['barcodes_unmapped'] = self.df_dict['barcodes'][-map_mask]
            self.df_dict['barcodes'] = self.df_dict['barcodes'][map_mask]
            del map_mask

        # count variants associated with the barcodes
        for bc, count in self.df_dict['barcodes'].iterrows():
            count = count['count']
            variant = self.barcode_map[bc]
            mutations = self.count_variant(variant, copies=count)
            if mutations is None: # variant has too many mutations
                self.filter_stats['max mutations'] += count
                self.filter_stats['total'] += count
                if self.verbose:
                    self.report_filtered_variant(variant, count)
            else:
                if mutations not in self.barcode_map.variants:
                    self.barcode_map.variants[mutations] = list()
                if bc not in self.barcode_map.variants[mutations]:
                    self.barcode_map.variants[mutations].append(bc)

        self.df_dict['variants'] = \
                pd.DataFrame.from_dict(self.df_dict['variants'], 
                                       orient="index", dtype="int32")
        if len(self.df_dict['variants']) == 0:
            raise EnrichError("Failed to count variants", self.name)
        self.df_dict['variants'].columns = ['count']
        self.df_dict['variants'].sort('count', ascending=False, inplace=True)

        logging.info("Counted %d variants (%d unique) [%s]" % \
                (self.df_dict['variants']['count'].sum(), len(self.df_dict['variants'].index), self.name))
        if self.aligner is not None:
            logging.info("Aligned %d variants [%s]" % (self.aligner.calls, self.name))
        self.report_filter_stats()


    def orphan_barcodes(self, mincount=0):
        """
        Returns a list of barcodes that are not found in the object's :py:class:`BarcodeMap` 
        but are present in the library more than *mincount* times.
        """
        orphans = [x in self.barcode_map for x in \
                   self.df_dict['barcodes'].index]
        return self.df_dict['barcodes'][orphans] \
                [self.df_dict['barcodes']['count'] >= mincount]


    def report_filtered_variant(self, variant, count):
        """
        Outputs a summary of the filtered variant to *handle*. Internal filter 
        names are converted to messages using the ``DataContainer._filter_messages`` 
        dictionary. Related to :py:meth:`SeqLib.report_filtered`.
        """
        logging.debug("Filtered variant (quantity=%d) (%s) [%s]\n%s" % \
                    (count, DataContainer._filter_messages['max mutations'], self.name, variant), file=handle)

