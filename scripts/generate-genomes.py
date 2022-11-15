#! /usr/bin/env python
import argparse
import sys
import random
import gzip

# some code taken from https://github.com/ctb/dbg-graph-null

x = ["A"] + ["G"] + ["C"] + ["T"]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--length', type=int, default=50000,
                        help="Simulated genome length")
    parser.add_argument('-n', '--number', type=int, default=1000,
                        help="number of sequences to generate")
    parser.add_argument('-o', '--output', required=True,
                        help="output sequences to this (gzipped) file")
    args = parser.parse_args()

    LENGTH = args.length
    assert LENGTH % 4 == 0, "length must be divisible by 4"

    with gzip.open(args.output, 'wt') as fp:
        sys.stderr.write("... 0")
        sys.stderr.flush()

        for i in range(1, args.number + 1):
            if i % 100 == 0:
                sys.stderr.write(u'\r\033[K')
                sys.stderr.write(f"... {i}")
                sys.stderr.flush()
            random.seed(i)

            dna = x * (LENGTH // 4)
            random.shuffle(dna)
            dna = "".join(dna)
        
            fp.write(f">genome_{i}\n{dna}\n")

    sys.stderr.write("\n")
    sys.stderr.write(F"Generated {args.number} sequences of length {LENGTH}\n")


if __name__ == '__main__':
    sys.exit(main())
