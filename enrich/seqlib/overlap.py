from __future__ import print_function
from sys import stderr
from variant import VariantSeqLib
from enrich_error import EnrichError
from fqread import read_fastq_multi, check_fastq, FQRead
import pandas as pd
import logging


class OverlapSeqLib(VariantSeqLib):
    """
    Class for count data from sequencing libraries with overlapping paired-end 
    reads for each variant. Creating a 
    :py:class:`~seqlib.overlap.OverlapSeqLib` requires a valid *config* object 
    with an ``'overlap'`` entry.

    Example config file for a :py:class:`OverlapSeqLib`:

    .. literalinclude:: config_examples/overlap.json

    :download:`Download this JSON file <config_examples/overlap.json>`

    The ``"fastq"`` config entry must contain two read files, with the keys 
    ``"forward"`` and ``"reverse"``. Information about how to combine these 
    reads is in the ``"overlap"`` config entry.

    The ``"overlap"`` config entry contains the following keys:

    * ``"forward start"`` --- position in the forward read where the \
        overlapping region begins
    * ``"reverse start"`` --- position in the reverse read where the \
        overlapping region begins (before being reverse-complemented)
    * ``"length"`` --- number of bases in the overlapping region
    * ``"max mismatches"`` --- maximum number of mismatches tolerated in the \
        overlapping region before discarding the read
    * ``"overlap only"`` --- whether to trim the merged read to contain only \
        the overlapping region (optional, defaults to ``False``)

    Here is a schematic of the case in the above JSON example::

        forward ---> 1   
                     CGACGCAAGGA
                       |||||||||
                       ACTCCTTGCGTCG
                                   1 <--- reverse

    Note that the merged sequence is identical to the wild type sequence given 
    in the JSON file.
    """
    def __init__(self, config):
        VariantSeqLib.__init__(self, config)
        try:
            self.forward = config['fastq']['forward']
            self.reverse = config['fastq']['reverse']

            self.fwd_start = int(config['overlap']['forward start'])
            self.rev_start = int(config['overlap']['reverse start'])
            self.overlap_length = int(config['overlap']['length'])
            self.trim = config['overlap']['overlap only']
            self.max_overlap_mismatches = int(config['overlap']
                                                    ['max mismatches'])

            if 'merge failure' in config['filters']:
                raise EnrichError("'merge failure' is not user-configurable", 
                                  self.name)
            self.set_filters(config['filters'], {'remove unresolvable' : False, 
                                      'min quality' : 0,
                                      'avg quality' : 0,
                                      'max mutations' : len(self.wt_dna),
                                      'chastity' : False,
                                      'merge failure' : True})
        except KeyError as key:
            raise EnrichError("Missing required config value {key}".format(key=key), 
                              self.name)
        except ValueError as value:
            raise EnrichError("Invalid parameter value {value}".format(value=value), self.name)

        try:
            check_fastq(self.forward)
            check_fastq(self.reverse)
        except IOError as fqerr:
            raise EnrichError("FASTQ file error: {error}".format(error=fqerr), self.name)


    def merge_reads(self, fwd, rev):
        """
        Combines the *fwd* and *rev* :py:class:`~fqread.FQRead` objects into a 
        single :py:class:`~fqread.FQRead` with the same header information as 
        *fwd*. Mismatches are resolved by taking the highest quality base. If 
        discrepant bases have the same quality value, this position is 
        unresolvable and an ``'X'`` is inserted. Quality values in the 
        resulting :py:class:`~fqread.FQRead` are the maximum quality for the 
        given base at that position. Returns ``None`` if the maximum number of 
        mismatches in the overlap region is exceded.
        """
        rev.revcomp()

        rev_extra_start = len(rev) - self.rev_start + 1
        fwd_end = self.fwd_start + self.overlap_length - 1
        merge = FQRead(header=fwd.header, 
                       sequence="A",
                       header2=fwd.header2,
                       quality="#",
                       qbase=fwd.qbase)
        merge.sequence = list(fwd.sequence[:fwd_end] + \
                                     rev.sequence[rev_extra_start:])
        merge.quality = fwd.quality[:fwd_end] + \
                                     rev.quality[rev_extra_start:]

        mismatches = 0
        for i in xrange(self.overlap_length):
            a = self.fwd_start - 1 + i
            b = len(rev) - self.rev_start - self.overlap_length + i + 1
            if fwd.sequence[a] == rev.sequence[b]:
                # take the highest quality value
                if rev.quality[b] > fwd.quality[a]:
                    merge.quality[a] = rev.quality[b]
            else:
                mismatches += 1
                if fwd.quality[a] == rev.quality[b]:
                    merge.sequence[a] = 'X' # unresolvable
                elif rev.quality[b] > fwd.quality[a]:
                    merge.sequence[a] = rev.sequence[b]
                    merge.quality[a] = rev.quality[b]
                else:
                    pass # overlap region already same as fwd
        if mismatches > self.max_overlap_mismatches:
            return None # merge failed

        merge.sequence = "".join(merge.sequence)
        if self.trim:
            merge.trim_length(self.overlap_length, self.fwd_start)
        return merge


    def calculate(self):
        """
        Reads the forward and reverse reads, merges them, performs 
        quality-based filtering, and counts the variants.
        """
        self.df_dict['variants'] = dict()

        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        logging.info("Counting variants [{name}]".format(name=self.name))
        for fwd, rev in read_fastq_multi([self.forward, self.reverse]):
            for key in filter_flags:
                filter_flags[key] = False

            # filter the read based on specified quality settings
            if self.filters['chastity']:
                if not fwd.is_chaste():
                    filter_flags['chastity'] = True
                    if self.report_filtered:
                        self.report_filtered_read(fwd, filter_flags)
                if not rev.is_chaste():
                    filter_flags['chastity'] = True
                    if self.report_filtered:
                        self.report_filtered_read(rev, filter_flags)
                if filter_flags['chastity']:
                    self.filter_stats['chastity'] += 1
                    self.filter_stats['total'] += 1
                    continue
            merge = self.merge_reads(fwd, rev)
            if merge is None: # merge failed
                self.filter_stats['merge failure'] += 1
                self.filter_stats['total'] += 1
                filter_flags['merge failure'] = True
                if self.report_filtered:
                    self.report_filtered_read(fwd, filter_flags)
                    self.report_filtered_read(rev, filter_flags)
            else:
                if self.filters['remove unresolvable']:
                    if 'X' in merge.sequence:
                        self.filter_stats['remove unresolvable'] += 1
                        filter_flags['remove unresolvable'] = True
                if self.filters['min quality'] > 0:
                    if merge.min_quality() < self.filters['min quality']:
                        self.filter_stats['min quality'] += 1
                        filter_flags['min quality'] = True
                if self.filters['avg quality'] > 0:
                    if merge.mean_quality() < self.filters['avg quality']:
                        self.filter_stats['avg quality'] += 1
                        filter_flags['avg quality'] = True
                if not any(filter_flags.values()): # passed quality filtering
                    mutations = self.count_variant(merge.sequence)
                    if mutations is None: # merge read has too many mutations
                        self.filter_stats['max mutations'] += 1
                        filter_flags['max mutations'] = True
                if any(filter_flags.values()):
                    self.filter_stats['total'] += 1
                    if self.report_filtered:
                        self.report_filtered_read(merge, filter_flags)

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

