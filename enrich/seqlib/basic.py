from variant import VariantSeqLib
from enrich_error import EnrichError
from fqread import read_fastq, check_fastq
import pandas as pd
import logging


class BasicSeqLib(VariantSeqLib):
    """
    Class for count data from sequencing libraries with a single read for 
    each variant. Creating a :py:class:`BasicSeqLib` requires a valid 
    *config* object, usually from a ``.json`` configuration file.
    """
    def __init__(self, config):
        VariantSeqLib.__init__(self, config)
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

            self.set_filters(config['filters'], {'min quality' : 0,
                                      'avg quality' : 0,
                                      'chastity' : False,
                                      'max mutations' : len(self.wt_dna)})
        except KeyError as key:
            raise EnrichError("missing required config value: {key}".format(key=key), self.name)

        try:
            check_fastq(self.reads)
        except IOError as fqerr:
            raise EnrichError("FASTQ file error: {error}".format(error=fqerr), self.name)


    def calculate(self):
        """
        Reads the forward or reverse FASTQ file (reverse reads are reverse-complemented),
        performs quality-based filtering, and counts the variants.
        """
        self.df_dict['variants'] = dict()

        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        logging.info("Counting variants [{name}]".format(name=self.name))
        for fq in read_fastq(self.reads):
            if self.revcomp_reads:
                fq.revcomp()

            for key in filter_flags:
                filter_flags[key] = False

            # filter the read based on specified quality settings
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
            if not any(filter_flags.values()): # passed quality filtering
                mutations = self.count_variant(fq.sequence)
                if mutations is None: # read has too many mutations
                    self.filter_stats['max mutations'] += 1
                    filter_flags['max mutations'] = True
            if any(filter_flags.values()):
                self.filter_stats['total'] += 1
                if self.report_filtered:
                    self.report_filtered_read(fq, filter_flags)

        self.df_dict['variants'] = \
                pd.DataFrame.from_dict(self.df_dict['variants'], 
                                       orient="index", dtype="int32")
        if len(self.df_dict['variants']) == 0:
            raise EnrichError("Failed to count variants", self.name)
        self.df_dict['variants'].columns = ['count']
        self.df_dict['variants'].sort('count', ascending=False, inplace=True)

        logging.info("Counted {n} variants ({u} unique) [{name}]".format(
                n=self.df_dict['variants']['count'].sum(), u=len(self.df_dict['variants'].index), name=self.name))
        if self.aligner is not None:
            logging.info("Aligned {n} variants [{name}]".format(n=self.aligner.calls, name=self.name))
            self.aligner_cache = None
        self.report_filter_stats()
