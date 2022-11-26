#! /usr/bin/env python
import sys
import screed
import random
import argparse

import sequtils                 # local import


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--error-rate', type=float, default=.01)
    parser.add_argument('-r', '--read-length', type=int, default=100,
                        help="Length of reads to generate")
    parser.add_argument('-C', '--coverage', type=int, default=50,
                        help="Targeted coverage level")
    parser.add_argument("-S", "--seed", dest="seed", help="Random seed", type=int,
                        default=1)

    parser.add_argument('genome')
    args = parser.parse_args()

    COVERAGE=args.coverage
    READLEN=args.read_length
    ERROR_RATE=args.error_rate

    random.seed(args.seed)                  # make this reproducible, please.

    records = list(screed.open(args.genome))
    
    len_genome = sum([ len(x.sequence) for x in records ])

    print('genome size:', len_genome, file=sys.stderr)
    print('coverage:', COVERAGE, file=sys.stderr)
    print('readlen:', READLEN, file=sys.stderr)
    print('error rate:', ERROR_RATE, file=sys.stderr)

    n_reads = int(len_genome*COVERAGE / float(READLEN))
    reads_mut = 0
    total_mut = 0

    print(f"Read in template genome of length {len_genome} from {args.genome}", file=sys.stderr)
    print(f"Generating {n_reads} reads of length {READLEN} for a target coverage of {COVERAGE} with a target error rate of 1 in {ERROR_RATE}", file=sys.stderr)

    for i, res in zip(range(n_reads), sequtils.generate_mutated_reads(records,
                                                                      READLEN,
                                                                      ERROR_RATE)):
        start, read, read_mutations = res

        if read_mutations:
            reads_mut += 1
            total_mut += read_mutations

        seq_name = f"read{i}"

        #print(f"@{seq_name} start={start},mutations={read_mutations}\n{read}\n+\n{'B'*READLEN}")
        print(f">{seq_name} start={start},mutations={read_mutations}\n{read}")

    print(f"{reads_mut} of {n_reads} reads mutated; {total_mut} total mutations", file=sys.stderr)


if __name__ == '__main__':
    sys.exit(main())
