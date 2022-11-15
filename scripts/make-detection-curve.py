#! /usr/bin/env python
import sys
import csv
import screed
import random
import argparse

import sourmash

import sequtils                 # local import


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('genome')
    parser.add_argument('-e', '--error-rate', type=float, default=.01)
    parser.add_argument('-r', '--read-length', type=int, default=100,
                        help="Length of reads to generate")
    parser.add_argument("-S", "--seed", dest="seed", help="Random seed", type=int,
                        default=1)
    parser.add_argument("-k", "--ksize", default=31, help="k-mer size")
    parser.add_argument("-o", "--output", required=True,
                        help="CSV output of detection curve")

    args = parser.parse_args()

    READLEN=args.read_length
    ERROR_RATE=args.error_rate
    NUM_FRACMINHASH = 5

    random.seed(args.seed)                  # make this reproducible, please.

    records = list(screed.open(args.genome))
    assert len(records) == 1
    record = records[0]
    
    genome = record.sequence
    len_genome = len(genome)

    total_mh = sourmash.MinHash(0, args.ksize, scaled=1)
    total_mh.add_sequence(genome)
    all_hashes = set(total_mh.hashes)

    # make NUM_FRACMINHASH minhashes each with different mmh seeds
    all_hashes_list = []
    scaled_mh_list = []
    for i in range(NUM_FRACMINHASH):
        smh = sourmash.MinHash(0, args.ksize, scaled=1000, seed=i + 42)
        all_hashes_i = smh.copy_and_clear()
        all_hashes_i.add_sequence(genome)
        scaled_mh_list.append(smh)
        all_hashes_list.append(all_hashes_i)
    
    print('genome size:', len_genome, file=sys.stderr)
    print('readlen:', READLEN, file=sys.stderr)
    print('error rate:', ERROR_RATE, file=sys.stderr)
    print('num k-mers:', len(total_mh))

    reads_mut = 0
    total_mut = 0

    print(f"Read in template genome {0} of length {1} from {2}".format(record["name"], len_genome, args.genome), file=sys.stderr)
    print(f"Generating reads of length {READLEN} with an error rate of 1 in {ERROR_RATE}", file=sys.stderr)

    it = sequtils.generate_mutated_reads(genome, READLEN, ERROR_RATE)
    it = iter(it)

    fp = open(args.output, 'w', newline="")
    csv_w = csv.writer(fp)
    headers = ['num_reads', 'coverage', 'n_detected', 'f_detected']
    for i in range(NUM_FRACMINHASH):
        headers.append(f"smash_count_{i}")
    csv_w.writerow(headers)
    csv_w.writerow([0, 0, 0, 0] + [0]*NUM_FRACMINHASH)

    n_reads = 0
    total_bp_in_reads = 0
    f01 = len(all_hashes) * 0.1

    remaining_hashes = set(all_hashes)
    while len(remaining_hashes) > f01:
        start, read, read_mutations = next(it)

        if read_mutations:
            reads_mut += 1
            total_mut += read_mutations

        n_reads += 1
        total_bp_in_reads += len(read)

        # first, track _all_ hashes for actual k-mer detection
        mh = total_mh.copy_and_clear()
        mh.add_sequence(read)
        remaining_hashes -= set(mh.hashes)

        n_detected = len(all_hashes) - len(remaining_hashes)
        f_detected = n_detected / len(all_hashes)
        coverage = total_bp_in_reads / len_genome

        # now, track sourmash detection & intersect with legit hashes:
        smash_detection = []
        for smh, all_hashes_i in zip(scaled_mh_list, all_hashes_list):
            smh.add_sequence(read)
            smh_hashes = set(smh.hashes)
            smh_hashes.intersection_update(all_hashes_i.hashes)
            smash_detection.append(len(smh_hashes))

        csv_w.writerow([n_reads, f"{coverage:.4f}", n_detected, f"{f_detected:.4f}"] + smash_detection)

        sys.stdout.write(u'\r\033[K')
        sys.stdout.write(f"...{n_reads} reads, {len(all_hashes)} missing k-mers, {total_bp_in_reads / len_genome:.2f} coverage")
        sys.stdout.flush()

    fp.close()


if __name__ == '__main__':
    sys.exit(main())
