#! /usr/bin/env python
"""
Given a query sig, and a target set of contigs,

extract the maximum extent of actual DNA sequence from the contigs that
would generate the hashes in the query sig. For example, for a single hash,
this would be all DNA sequence between the neighboring hashes, but NOT
including the k-mers that generated those neighboring hashes. If multiple
hashes in a row matched, it would include all sequence up to the surrounding
(non-matching) k-mers.

The motivation is to extract the maximum potentially alignable sequence
when two genomes share some hashes, but might be useful for other purposes,
too.

Warning: not yet well tested! Smoketests seem to pass :shrug:
"""
import sys
import argparse
import screed
import sourmash
from collections import defaultdict, Counter
from interval import interval
import gzip

# QUESTIONS: do hashes tend to be in runs?
# TODO: add tests...
# * query = genome => all same

def extract_runs(hash_positions, len_sequence, *, min_runcount=1):
    ivals = []

    # special case first
    last_was_in_query = False
    runcount = 0
    if hash_positions[0][1]:
        run_start_pos = 0
        last_was_in_query = True
        runcount = 1

    i = 1
    while i < len(hash_positions):
        pos, in_query, *rest = hash_positions[i]

        if in_query:
            if last_was_in_query:
                # last was in query, this is in query
                # CONTINUE RUN
                runcount += 1
            else: # start run at last position + 1
                run_start_pos = hash_positions[i-1][0] + 1
                last_was_in_query = True
                runcount = 1
        else:
            if last_was_in_query: # end run at this position - 1
                if runcount >= min_runcount:
                    ivals.append(interval([run_start_pos, pos - 1]))
                last_was_in_query = False
                run_start_pos = None
            else:
                # not in query, last was not in query
                # NO RUN TO START/CONTINUE
                pass

        i += 1

    # special case end
    if last_was_in_query:
        if runcount >= min_runcount:
            ivals.append(interval[run_start_pos, len_sequence - 1])

    return ivals


def main():
    p = argparse.ArgumentParser()
    p.add_argument('contigs')
    p.add_argument('sigfile')
    p.add_argument('-o', '--output-contigs')
    args = p.parse_args()

    idx = sourmash.load_file_as_index(args.sigfile)
    idx = idx.select(ksize=31, moltype='DNA')
    sigs = list(idx.signatures())
    assert len(sigs) == 1
    query_sig = sigs[0]
    query_hashes = set(query_sig.minhash.hashes)

    filtered_seqs = []
    for record in screed.open(args.contigs):
        # calculate basic FracMinHash
        background_mh = query_sig.minhash.copy_and_clear()
        background_mh.add_sequence(record.sequence, force=True)
        background_hashes = set(background_mh.hashes)
        
        # does this sequence have any relevant hashes in it?
        if not background_hashes.intersection(query_hashes):
            # no, ignore!
            continue
        else:
            # go through _all_ kmers/hashes, not just fracminhash subset
            hash_iter = background_mh.seq_to_hashes(record.sequence,
                                                    force=True,
                                                    bad_kmers_as_zeroes=True)

            hash_positions = []
            for pos, hashval in enumerate(hash_iter):
                # first, filter on whether it belongs in sketch
                if hashval and hashval in background_hashes:
                    in_query = False

                    # now, figure out if it's one of our query hashes
                    if hashval in query_hashes:
                        in_query = True

                    # record!
                    hash_positions.append((pos, in_query, hashval))

            # now, convert to intervals
            runs = extract_runs(hash_positions, len(record.sequence))

            # now... combine?
            #print(runs)

            total_run_seq = 0
            for ival in runs:
                start, end = map(int, ival[0])
                total_run_seq += end - start + 1
                filtered_seqs.append((record.name, record.sequence[start:end + 1]))

            print(record.name, len(record.sequence), total_run_seq)

    # double check?
    check_mh = query_sig.minhash.copy_and_clear()
    for (name, seq) in filtered_seqs:
        check_mh.add_sequence(seq)

    # note: there may be hashes in query that are NOT IN this genome!
    assert check_mh.contained_by(query_sig.minhash) == 1.0

    # output?
    if args.output_contigs:
        xopen = open
        if args.output_contigs.endswith('.gz'):
            xopen = gzip.open
        with xopen(args.output_contigs, 'wt') as fp:
            for n, (name, seq) in enumerate(filtered_seqs):
                fp.write(f">{name}.{n}\n{seq}\n")


if __name__ == '__main__':
    sys.exit(main())
