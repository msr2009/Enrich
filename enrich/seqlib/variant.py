from __future__ import print_function
import re
from sys import stdout, stderr
from enrich_error import EnrichError
from aligner import Aligner
from seqlib import SeqLib
import pandas as pd


# Variant string for counting wild type sequences
WILD_TYPE_VARIANT = "_wt"

# Standard codon table for translating wild type and variant DNA sequences
codon_table = {
        'TTT':'F', 'TCT':'S', 'TAT':'Y', 'TGT':'C',
        'TTC':'F', 'TCC':'S', 'TAC':'Y', 'TGC':'C',
        'TTA':'L', 'TCA':'S', 'TAA':'*', 'TGA':'*',
        'TTG':'L', 'TCG':'S', 'TAG':'*', 'TGG':'W',
        'CTT':'L', 'CCT':'P', 'CAT':'H', 'CGT':'R',
        'CTC':'L', 'CCC':'P', 'CAC':'H', 'CGC':'R',
        'CTA':'L', 'CCA':'P', 'CAA':'Q', 'CGA':'R',
        'CTG':'L', 'CCG':'P', 'CAG':'Q', 'CGG':'R',
        'ATT':'I', 'ACT':'T', 'AAT':'N', 'AGT':'S',
        'ATC':'I', 'ACC':'T', 'AAC':'N', 'AGC':'S',
        'ATA':'I', 'ACA':'T', 'AAA':'K', 'AGA':'R',
        'ATG':'M', 'ACG':'T', 'AAG':'K', 'AGG':'R',
        'GTT':'V', 'GCT':'A', 'GAT':'D', 'GGT':'G',
        'GTC':'V', 'GCC':'A', 'GAC':'D', 'GGC':'G',
        'GTA':'V', 'GCA':'A', 'GAA':'E', 'GGA':'G',
        'GTG':'V', 'GCG':'A', 'GAG':'E', 'GGG':'G'
    }


# Convert between single- and three-letter amino acid codes
aa_codes = {
        'Ala' : 'A', 'A' : 'Ala',
        'Arg' : 'R', 'R' : 'Arg',
        'Asn' : 'N', 'N' : 'Asn',
        'Asp' : 'D', 'D' : 'Asp',
        'Cys' : 'C', 'C' : 'Cys',
        'Glu' : 'E', 'E' : 'Glu',
        'Gln' : 'Q', 'Q' : 'Gln',
        'Gly' : 'G', 'G' : 'Gly',
        'His' : 'H', 'H' : 'His',
        'Ile' : 'I', 'I' : 'Ile',
        'Leu' : 'L', 'L' : 'Leu',
        'Lys' : 'K', 'K' : 'Lys',
        'Met' : 'M', 'M' : 'Met',
        'Phe' : 'F', 'F' : 'Phe',
        'Pro' : 'P', 'P' : 'Pro',
        'Ser' : 'S', 'S' : 'Ser',
        'Thr' : 'T', 'T' : 'Thr',
        'Trp' : 'W', 'W' : 'Trp',
        'Tyr' : 'Y', 'Y' : 'Tyr',
        'Val' : 'V', 'V' : 'Val',
        'Ter' : '*', '*' : 'Ter',
        '?' : '???'
}


def has_indel(variant):
    """
    Helper function that returns ``True`` if the HGVS string *variant* 
    contains an indel mutation.
    """
    return any(x in variant for x in ("ins", "del", "dup"))


class VariantSeqLib(SeqLib):
    """
    Abstract :py:class:`SeqLib` class for for Enrich libraries containing variants. Implements core functionality for assessing variants, either coding
    or noncoding. Subclasess must evaluate the variant DNA sequences that are being counted.
    """
    def __init__(self, config, parent=True):
        if parent:
            SeqLib.__init__(self, config)
        self.wt_dna = None
        self.wt_protein = None
        self.aligner = None

        try:
            self.set_wt(config['wild type']['sequence'], 
                        coding=config['wild type']['coding'])
            if 'align variants' in config:
                if config['align variants']:
                    self.aligner = Aligner()
                else:
                    self.aligner = None
            else:
                self.aligner = None

        except KeyError as key:
            raise EnrichError("Missing required config value '{key}'".format(key), 
                              self.name)

        if 'reference offset' in config['wild type']:
            try:
                self.reference_offset = int(config['wild type']
                                                  ['reference offset'])
            except ValueError:
                raise EnrichError("Invalid reference offset value", self.name)
        else:
            self.reference_offset = 0

        self.df_dict['variants'] = None


    def is_coding(self):
        return self.wt_protein is not None


    def set_wt(self, sequence, coding=True):
        """
        Set the wild type DNA *sequence*. The *sequence* is translated if *coding* 
        is ``True``. The *sequence* may only contain ``ACGT``, but may 
        contain whitespace (which will be removed). If *coding*, *sequence* must be in-frame.
        """
        sequence = "".join(sequence.split()) # remove whitespace

        if not re.match("^[ACGTacgt]+$", sequence):
            raise EnrichError("WT DNA sequence contains unexpected "
                              "characters", self.name)
        if len(sequence) % 3 != 0 and coding:
            raise EnrichError("WT DNA sequence contains incomplete codons", 
                              self.name)
        
        self.wt_dna = sequence.upper()
        if coding:
            self.wt_protein = ""
            for i in xrange(0, len(self.wt_dna), 3):
                self.wt_protein += codon_table[self.wt_dna[i:i + 3]]
        else:
            self.wt_protein = None


    def align_variant(self, variant_dna):
        """
        Use the local :py:class:`~seqlib.aligner.Aligner` instance to align the *variant_dna* to the 
        wild type sequence. Returns a list of HGVS variant strings.

        .. warning:: Using the :py:class:`~seqlib.aligner.Aligner` dramatically increases runtime.
        """
        mutations = list()
        traceback = self.aligner.align(self.wt_dna, variant_dna)
        for x, y, cat, length in traceback:
            if cat == "match":
                continue
            elif cat == "mismatch":
                mut = "{pre}>{post}".format(pre=self.wt_dna[x], post=variant_dna[y])
            elif cat == "insertion":
                if y > length:
                    dup = variant_dna[y:y + length]
                    if dup == variant_dna[y - length:y]:
                        mut = "dup{seq}".format(seq=dup)
                    else:
                        mut = "_{pos}ins{seq}".format(post=x + 2, seq=dup)
                else:                                    
                    mut = "_{pos}ins{seq}".format(pos=x + 2, seq=variant_dna[y:y + length])
            elif cat == "deletion":
                mut = "_{pos}del".format(pos=x + length)
            mutations.append((x, mut))
        return mutations


    def count_variant(self, variant_dna, copies=1, include_indels=True):
        """
        Identifies mutations and counts the *variant_dna* sequence.
        The algorithm attempts to call variants by comparing base-by-base.
        If the *variant_dna* and wild type DNA are different lengths, or if there
        are an excess of mismatches (indicating a possible indel), local
        alignment is performed using :py:meth:`align_variant` if this option 
        has been selected in the configuration.

        Each variant is stored as a tab-delimited string of mutations in HGVS 
        format. Returns a list of HGSV variant strings. Returns an empty list 
        if the variant is wild type. Returns None if the variant was discarded
        due to excess mismatches.
        """
        if not re.match("^[ACGTNXacgtnx]+$", variant_dna):
            raise EnrichError("Variant DNA sequence contains unexpected "
                              "characters", self.name)

        variant_dna = variant_dna.upper()

        if len(variant_dna) != len(self.wt_dna):
            if self.aligner is not None:
                mutations = self.align_variant(variant_dna)
            else:
                return None
        else:
            mutations = list()
            for i in xrange(len(variant_dna)):
                if variant_dna[i] != self.wt_dna[i]:
                    mutations.append((i, "{pre}>{post}".format(pre=self.wt_dna[i], post=variant_dna[i])))
                    if len(mutations) > self.filters['max mutations']:
                        if self.aligner is not None:
                            mutations = self.align_variant(variant_dna)
                            if len(mutations) > self.filters['max mutations']:
                                # too many mutations post-alignment
                                return None
                            else:
                                # stop looping over this variant
                                break
                        else:
                            # too many mutations and not using aligner
                            return None

        mutation_strings = list()
        if self.is_coding():
            variant_protein = ""
            for i in xrange(0, len(variant_dna), 3):
                try:
                    variant_protein += codon_table[variant_dna[i:i + 3]]
                except KeyError: # garbage codon due to indel
                    variant_protein += '?'

            for pos, change in mutations:
                ref_dna_pos = pos + self.reference_offset + 1
                ref_pro_pos = (pos + self.reference_offset) / 3 + 1
                mut = "c.{pos}{change}".format(pos=ref_dna_pos, change=change)
                if has_indel(change):
                    mut += " (p.{pre}{pos}fs)".format(pre=aa_codes[self.wt_protein[pos / 3]], pos=ref_pro_pos)
                elif variant_protein[pos / 3] == self.wt_protein[pos / 3]:
                    mut += " (p.=)"
                else:
                    mut += " (p.{pre}{pos}{post})".format(pre=aa_codes[self.wt_protein[pos / 3]], pos=ref_pro_pos,
                             post=aa_codes[variant_protein[pos / 3]])
                mutation_strings.append(mut)
        else:
            for pos, change in mutations:
                ref_dna_pos = pos + self.reference_offset + 1
                mut = "n.{pos}{change}".format(pos=ref_dna_pos, change=change)
                mutation_strings.append(mut)

        if len(mutation_strings) > 0:
            variant_string = ', '.join(mutation_strings)
        else:
            variant_string = WILD_TYPE_VARIANT
        try:
            self.df_dict['variants'][variant_string] += copies
        except KeyError:
            self.df_dict['variants'][variant_string] = copies
        return variant_string


    def count_mutations(self, include_indels=False):
        """
        Count the individual mutations in all variants. If *include_indels* is ``False``, all mutations in a variant that contains 
        an insertion/deletion/duplication will not be counted. For coding sequences, amino acid substitutions are counted
        independently of the corresponding nucleotide change.
        """
        # restore the counts if they were saved to disk
        if self.df_dict['variants'] is None:
            self.load_counts(keys=['variants'])

        # create new dictionaries
        self.df_dict['mutations_nt'] = dict()
        if self.is_coding():
            self.df_dict['mutations_aa'] = dict()

        if not include_indels:
            mask = self.df_dict['variants'].index.map(has_indel)
            variant_data = self.df_dict['variants'][np.invert(mask)]
            del mask
        else:
            variant_data = self.df_dict['variants']
        if self.is_coding():
            for variant, count in variant_data.iterrows():
                count = count['count'] # get the element from the Series
                mutations = variant.split(", ")
                # get just the nucleotide changes
                for m in mutations:
                    m = m.split(" (")[0]
                    try:
                        self.df_dict['mutations_nt'][m] += count
                    except KeyError:
                        self.df_dict['mutations_nt'][m] = count
                # get the amino acid changes
                aa_changes = re.findall("p\.[A-Z][a-z][a-z]\d+[A-Z][a-z][a-z]", variant)
                for a in aa_changes:
                    try:
                        self.df_dict['mutations_aa'][a] += count
                    except KeyError:
                        self.df_dict['mutations_aa'][a] = count
        else:
            for variant, count in variant_data.iterrows():
                count = count['count'] # get the element from the Series
                mutations = variant.split(", ")
                for m in mutations:
                    try:
                        self.df_dict['mutations_nt'][m] += count
                    except KeyError:
                        self.df_dict['mutations_nt'][m] = count

        self.df_dict['mutations_nt'] = \
                pd.DataFrame.from_dict(self.df_dict['mutations_nt'], 
                                       orient="index", dtype="int32")
        if self.is_coding():
            self.df_dict['mutations_aa'] = \
                    pd.DataFrame.from_dict(self.df_dict['mutations_aa'], 
                                           orient="index", dtype="int32")

