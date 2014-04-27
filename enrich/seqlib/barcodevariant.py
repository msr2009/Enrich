from __future__ import print_function
import re
import logging
import gzip
import bz2
import os.path
from variant import VariantSeqLib
from barcode import BarcodeSeqLib
from seqlib import SeqLib
from enrich_error import EnrichError
from fqread import read_fastq, check_fastq
import pandas as pd


# variant string for a filtered variant
FILTERED_VARIANT = "~filtered"


class BarcodeMap(dict):
    """
    Dictionary-derived class for storing the relationship between barcodes 
    and variants. Requires the path to a *mapfile*, containing lines in the 
    format ``'barcode<tab>variant'`` for each barcode expected in the library. 
    This file can be plain text or compressed (``.bz2`` or ``.gz``).

    Barcodes must only contain the characters ``ACGT`` and variants must only 
    contain the characters ``ACGTN`` (lowercase characters are also accepted). 
    """
    def __init__(self, mapfile):
        self.name = "barcodemap_{fname}".format(fname=os.path.basename(mapfile))
        self.filename = mapfile
        self.variants = dict()
        self.bc_variant_strings = dict()

        try:
            ext = os.path.splitext(mapfile)[-1].lower()
            if ext in (".bz2"):
                handle = bz2.BZ2File(mapfile, "rU")
            elif ext in (".gz"):
                handle = gzip.GzipFile(mapfile, "rU")
            else:
                handle = open(mapfile, "rU")
        except IOError:
            raise EnrichError("Could not open barcode map file '{fname}'".format(fname=mapfile), self.name)

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
                    raise EnrichError("Barcode '{bc}' assigned to multiple unique variants".format(bc=barcode), self.name)
            else:
                self[barcode] = variant
        handle.close()


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
            raise EnrichError("Missing required config value {key}".format(key=key), 
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

        logging.info("Converting barcodes to variants [{name}]".format(name=self.name))
        if self.filter_unmapped:
            map_mask = self.df_dict['barcodes'].index.isin(self.barcode_map)
            self.df_dict['barcodes_unmapped'] = self.df_dict['barcodes'][-map_mask]
            self.df_dict['barcodes'] = self.df_dict['barcodes'][map_mask]
            del map_mask
            logging.info("Writing counts for {n} unique unmapped barcodes to disk [{name}]".format(n=len(self.df_dict['barcodes_unmapped']), name=self.name))
            self.dump_data(keys=['barcodes_unmapped']) # save memory

        # count variants associated with the barcodes
        for bc, count in self.df_dict['barcodes'].iterrows():
            count = count['count']
            variant = self.barcode_map[bc]
            mutations = self.count_variant(variant, copies=count)
            if mutations is None: # variant has too many mutations
                self.filter_stats['max mutations'] += count
                self.filter_stats['total'] += count
                if self.report_filtered:
                    self.report_filtered_variant(variant, count)
                if bc not in self.barcode_map.bc_variant_strings:
                    self.barcode_map.bc_variant_strings[bc] = FILTERED_VARIANT
            else:
                if mutations not in self.barcode_map.variants:
                    self.barcode_map.variants[mutations] = set()
                self.barcode_map.variants[mutations].update([bc])
                self.barcode_map.bc_variant_strings[bc] = mutations


        self.df_dict['variants'] = \
                pd.DataFrame.from_dict(self.df_dict['variants'], 
                                       orient="index", dtype="int32")
        if len(self.df_dict['variants']) == 0:
            raise EnrichError("Failed to count variants", self.name)
        self.df_dict['variants'].columns = ['count']
        self.df_dict['variants'].sort('count', ascending=False, inplace=True)

        logging.info("Retained counts for {n} variants ({u} unique) [{name}]".format(
                n=self.df_dict['variants']['count'].sum(), u=len(self.df_dict['variants'].index), name=self.name))
        if self.aligner is not None:
            logging.info("Aligned {n} variants [{name}]".format(n=self.aligner.calls, name=self.name))
            self.aligner_cache = None
        self.report_filter_stats()


    def report_filtered_variant(self, variant, count):
        """
        Outputs a summary of the filtered variant to *handle*. Internal filter 
        names are converted to messages using the ``DataContainer._filter_messages`` 
        dictionary. Related to :py:meth:`SeqLib.report_filtered`.
        """
        logging.debug("Filtered variant (quantity={n}) ({messages}) [{name}]\n{read!s}".format(
                    n=count, messages=DataContainer._filter_messages['max mutations'], name=self.name, read=variant), file=handle)

