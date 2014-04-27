from __future__ import print_function
import re
import logging
from seqlib import SeqLib
from enrich_error import EnrichError
from fqread import read_fastq, check_fastq
import pandas as pd

# debugging
from sys import stdout, stderr


class BarcodeSeqLib(SeqLib):
    """
    Class for count data from barcoded sequencing libraries. Designed for 
    barcode-only quantification or as a parent class for 
    :py:class:`~seqlib.barcodevariant.BarcodeVariantSeqLib`. Creating a 
    :py:class:`~seqlib.barcode.BarcodeSeqLib` requires a valid *config* 
    object with a ``"barcodes"`` entry (this entry can be empty).

    Example config file for a :py:class:`~seqlib.barcode.BarcodeSeqLib`:

    .. literalinclude:: config_examples/barcode.json

    :download:`Download this JSON file <config_examples/barcode.json>`

    The ``"fastq"`` config entry can contain one read file, with the key 
    ``"forward"`` or ``"reverse"``. If the read file is ``"reverse"``, all 
    barcodes will be reverse-complemented before being counted. The 
    ``"fastq"`` entry can contain optional values ``"start"`` and 
    ``"length"``, which will be used to trim the barcodes before counting. 
    Bases are counted starting at 1.

    The ``"min count"`` entry in ``"barcodes"`` is used to filter out 
    low-abundance barcodes (likely the result of technical artifacts). 
    Barcodes that occur less than ``"min count"`` times are removed from the 
    dataset. Setting the ``"min count"`` option appropriately can 
    dramatically improve execution time and reduce memory usage.

    .. note:: The :py:class:`~seqlib.barcodevariant.BarcodeVariantSeqLib` \
    class implements an alternative method for removing artifactual barcodes \
    that may be more appropriate for users of that module.
    """
    def __init__(self, config, barcodevariant=False):
        self.barcodevariant = barcodevariant
        if not self.barcodevariant:
            SeqLib.__init__(self, config)
        try:
            if 'forward' in config['fastq'] and 'reverse' in config['fastq']:
                raise EnrichError("Multiple FASTQ files specified", self.name)
            elif 'forward' in config['fastq']:
                self.reads = config['fastq']['forward']
                self.revcomp_reads = False
            elif 'reverse' in config['fastq']:
                self.reads = config['fastq']['reverse']
                self.revcomp_reads = True
            else:
                raise KeyError("'forward' or 'reverse'")

            if 'start' in config['fastq']:
                self.bc_start = config['fastq']['start']
            else:
                self.bc_start = 1
            if 'length' in config['fastq']:
                self.bc_length = config['fastq']['length']
            else:
                self.bc_length = 2147483647 # longer than any read... for now

            if 'min count' in config['barcodes']:
                self.min_count = config['barcodes']['min count']
            else:
                self.min_count = 0

            self.set_filters(config['filters'], {'min quality' : 0,
                                      'avg quality' : 0,
                                      'chastity' : False})
        except KeyError as key:
            raise EnrichError("Missing required config value {key}".format(key=key), self.name)

        try:
            check_fastq(self.reads)
        except IOError as fqerr:
            raise EnrichError("FASTQ file error: {error}".format(error=fqerr), self.name)

        self.df_dict['barcodes'] = None
        if self.min_count > 0:
            self.df_dict['barcodes_low_abundance'] = None


    def calculate(self):
        """
        Reads the forward or reverse FASTQ file (reverse reads are 
        reverse-complemented), performs quality-based filtering, and counts 
        the barcodes.
        """
        self.df_dict['barcodes'] = dict()

        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        # count all the barcodes
        logging.info("Counting barcodes [{name}]".format(name=self.name))
        for fq in read_fastq(self.reads):
            fq.trim_length(self.bc_length, start=self.bc_start)
            if self.revcomp_reads:
                fq.revcomp()

            for key in filter_flags:
                filter_flags[key] = False

            # filter the barcode based on specified quality settings
            if self.filters['chastity']:
                if not fq.is_chaste():
                    self.filter_stats['chastity'] += 1
                    filter_flags['chastity'] = True
            if self.filters['min quality'] > 0:
                if fq.min_quality() < self.filters['min quality']:
                    self.filter_stats['min quality'] += 1
                    filter_flags['min quality'] = True
            if self.filters['avg quality'] > 0:
                if fq.mean_quality() < self.filters['avg quality']:
                    self.filter_stats['avg quality'] += 1
                    filter_flags['avg quality'] = True
            if any(filter_flags.values()): # failed quality filtering
                self.filter_stats['total'] += 1
                if self.report_filtered:
                    self.report_filtered_read(fq, filter_flags)
            else: # passed quality filtering
                try:
                    self.df_dict['barcodes'][fq.sequence.upper()] += 1
                except KeyError:
                    self.df_dict['barcodes'][fq.sequence.upper()] = 1

        self.df_dict['barcodes'] = \
                pd.DataFrame.from_dict(self.df_dict['barcodes'], 
                                       orient="index", dtype="int32")
        if len(self.df_dict['barcodes']) == 0:
            raise EnrichError("Failed to count barcodes", self.name)
        self.df_dict['barcodes'].columns = ['count']
        self.df_dict['barcodes'].sort('count', ascending=False, inplace=True)
        if 'barcodes_low_abundance' in self.df_dict: # min count is set
            self.df_dict['barcodes_low_abundance'] = self.df_dict['barcodes'][self.df_dict['barcodes']['count'] < self.min_count]
            logging.info("Writing counts for {n} unique low-abundance barcodes to disk [{name}]".format(n=len(self.df_dict['barcodes_low_abundance']), name=self.name))
            self.dump_data(keys=['barcodes_low_abundance'])
            self.df_dict['barcodes'] = self.df_dict['barcodes'][self.df_dict['barcodes']['count'] >= self.min_count]

        logging.info("Retained counts for {n} barcodes ({u} unique) [{name}]".format(
                n=self.df_dict['barcodes']['count'].sum(), u=len(self.df_dict['barcodes'].index), name=self.name))
        if not self.barcodevariant:
            self.report_filter_stats()

