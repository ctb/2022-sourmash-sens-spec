#! /usr/bin/env python
import sourmash
import sys
import argparse
import screed

from sourmash.command_sketch import _signatures_for_sketch_factory


def consume(sig, filename):
    print(f"consuming '{filename}'")
    for record in screed.open(filename):
        sig.add_sequence(record.sequence, force=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query_genome')
    p.add_argument('against_genomes', nargs='+')
    p.add_argument('-k', '--ksize', default=31, type=int)
    args = p.parse_args()

    query_filename = args.query_genome
    against_filenames = set(args.against_genomes)
    if query_filename in against_filenames:
        print(f"Removing query '{query_filename}' from against list.")
        against_filenames.remove(query_filename)

    param_str = f'k={args.ksize},scaled=1'
    print("using:", param_str)

    factory = _signatures_for_sketch_factory([param_str], 'dna')
    sigs = factory()
    assert len(sigs) == 1
    query_sig = sigs[0]

    assert query_sig.minhash.ksize == args.ksize
    assert query_sig.minhash.scaled == 1

    consume(query_sig, query_filename)
    print(len(query_sig.minhash))

    query_set = set(query_sig.minhash.hashes)

    against_sigs = []
    for filename in against_filenames:
        sig = factory()[0]
        consume(sig, filename)

        query_set -= set(sig.minhash.hashes)
        print(len(sig.minhash), len(query_set))

    print(len(query_set))


if __name__ == '__main__':
    sys.exit(main())
