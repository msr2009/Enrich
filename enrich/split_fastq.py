from __future__ import print_function
from sys import stderr
import argparse
import os.path
from fqread import read_fastq_multi


def split_fastq(outdir, sequences, index, forward, reverse, max_mismatches):
    """
    """
    if index is None:
        print("Error: no index file specified for split_fastq", file=stderr)
        return

    if len(sequences) == 0:
        print("Error: no index sequences provided", file=stderr)
        return

    # build an iterator to process the files in parallel
    fq_handles = dict() # output file handles and index read sequences
    if forward is not None and reverse is not None:
        fq_iterator = read_fastq_multi([index, forward, reverse], 
                                       match_lengths=True)
        for s in sequences:
            name, ext = os.path.splitext(os.path.basename(index))
            index_name = "{name}_{seq}{ext}".format(name=name, seq=s, ext=ext)
            index_name = os.path.join(outdir, index_name)

            name, ext = os.path.splitext(os.path.basename(forward))
            forward_name = "{name}_{seq}{ext}".format(name=name, seq=s, 
                                                      ext=ext)
            forward_name = os.path.join(outdir, forward_name)

            name, ext = os.path.splitext(os.path.basename(reverse))
            reverse_name = "{name}_{seq}{ext}".format(name=name, seq=s, 
                                                      ext=ext)
            reverse_name = os.path.join(outdir, reverse_name)

            fq_handles[s] = \
                (open(index_name, "w"), open(forward_name, "w"),
                 open(reverse_name, "w"))
    elif forward is not None:
        fq_iterator = read_fastq_multi([index, forward], match_lengths=True)
        for s in sequences:
            name, ext = os.path.splitext(os.path.basename(index))
            index_name = "{name}_{seq}{ext}".format(name=name, seq=s, ext=ext)
            index_name = os.path.join(outdir, index_name)

            name, ext = os.path.splitext(os.path.basename(forward))
            forward_name = "{name}_{seq}{ext}".format(name=name, seq=s, 
                                                      ext=ext)
            forward_name = os.path.join(outdir, forward_name)

            fq_handles[s] = \
                (open(index_name, "w"), open(forward_name, "w"))
    elif reverse is not None:
        fq_iterator = read_fastq_multi([index, reverse], match_lengths=True)
        for s in sequences:
            name, ext = os.path.splitext(os.path.basename(index))
            index_name = "{name}_{seq}{ext}".format(name=name, seq=s, ext=ext)
            index_name = os.path.join(outdir, index_name)

            name, ext = os.path.splitext(os.path.basename(reverse))
            reverse_name = "{name}_{seq}{ext}".format(name=name, seq=s, 
                                                      ext=ext)
            reverse_name = os.path.join(outdir, reverse_name)

            fq_handles[s] = \
                (open(index_name, "w"), open(reverse_name, "w"))
    else:
        print("Error: no forward or reverse files specified for split_fastq",
              file=stderr)
        return

    for t in fq_iterator:
        if t is None:
            print("Warning: FASTQ files are not the same length", file=stderr)
            break
        index_sequence = t[0].sequence

        match = None
        for s in sequences:
            mismatches = 0
            for i in xrange(len(s)):
                if index_sequence[i] != s[i]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break
            if mismatches <= max_mismatches:
                match = s
                break

        if match:
            for i in xrange(len(t)):
                print(t[i], file=fq_handles[match][i])

    # close all the files
    for handle_tuple in fq_handles.values():
        for h in handle_tuple:
            h.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create new FASTQ files for "
                                     "reads with the given index sequences.")
    parser.add_argument("sequences", metavar="SEQ", nargs="+",
                        help="index sequence to match")
    parser.add_argument("-f", "--forward", metavar="FQ", 
                        help="forward read FASTQ file")
    parser.add_argument("-r", "--reverse", metavar="FQ", 
                        help="reverse read FASTQ file")
    parser.add_argument("-i", "--index", metavar="FQ", 
                        help="index read FASTQ file")
    parser.add_argument("-o", "--output", default=".", metavar="DIR",
                        help="output directory")
    parser.add_argument("-m", "--mismatches", metavar="N",
                        help="index read mismatch threshold",
                        type=int, default=0)

    args = parser.parse_args()

    split_fastq(args.output, args.sequences, args.index, args.forward, 
                args.reverse, args.mismatches)
    