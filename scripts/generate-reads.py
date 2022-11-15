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
    assert len(records) == 1
    record = records[0]
    
    genome = record.sequence
    len_genome = len(genome)

    print('genome size:', len_genome, file=sys.stderr)
    print('coverage:', COVERAGE, file=sys.stderr)
    print('readlen:', READLEN, file=sys.stderr)
    print('error rate:', ERROR_RATE, file=sys.stderr)

    n_reads = int(len_genome*COVERAGE / float(READLEN))
    reads_mut = 0
    total_mut = 0

    nucl = ['a', 'c', 'g', 't']

    print("Read in template genome {0} of length {1} from {2}".format(record["name"], len_genome, args.genome), file=sys.stderr)
    print("Generating {0} reads of length {1} for a target coverage of {2} with a target error rate of 1 in {3}".format(n_reads, READLEN, COVERAGE, ERROR_RATE), file=sys.stderr)

    for i, res in zip(range(n_reads), sequtils.generate_mutated_reads(genome,
                                                                      READLEN,
                                                                      ERROR_RATE)):
        start, read, read_mutations = res

        if read_mutations:
            reads_mut += 1
            total_mut += read_mutations

        seq_name = f"read{i}"

        print('>{0} start={1},mutations={2}\n{3}'.format(seq_name, start, read_mutations, read), file=sys.stderr)

    print("%d of %d reads mutated; %d total mutations" % \
          (reads_mut, n_reads, total_mut), file=sys.stderr)


if __name__ == '__main__':
    sys.exit(main())
