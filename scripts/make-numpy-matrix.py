#! /usr/bin/env python
import sys
import numpy as np
import argparse
import os
import screed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('genome_files', nargs='+')
    p.add_argument('-r', '--readlen', default=100, help='readlen')
    p.add_argument('-C', '--coverage', default=1, help='coverage')
    p.add_argument('-o', '--output-prefix', required=True)
    args = p.parse_args()

    genome_files = []
    genome_idx = {}
    genome_names = {}
    for n, orig_g in enumerate(args.genome_files):
        assert orig_g.endswith('.fna.gz')
        g = orig_g[:-7]
        assert not g.endswith('.')
        genome_files.append(g)
        genome_idx[g] = n
        genome_names[g] = orig_g

    readcounts = {}
    for g in genome_files:
        readfile = f"reads/{g}.reads.l{args.readlen}.C{args.coverage}.fa.gz"
        with screed.open(readfile) as it:
            readcounts[g] = len(list(it))

    mat = np.ones((len(genome_files), len(genome_files)))

    for g1 in genome_files:
        g1_idx = genome_idx[g1]
        for g2 in genome_files:
            g2_idx = genome_idx[g2]
            countfile = f"counts/{g1}.x.{g2}.l{args.readlen}.C{args.coverage}.count"
            assert os.path.exists(countfile), countfile
            with open(countfile, "rt") as fp:
                mapcount = int(fp.read())
            readcount = readcounts[g1]

            f = mapcount / readcount

            mat[g2_idx, g1_idx] = f

    print(f"saving to {args.output_prefix} and {args.output_prefix}.labels.txt")
    with open(args.output_prefix, "wb") as fp:
        np.save(fp, mat)

    with open(f"{args.output_prefix}.labels.txt", "wt") as fp:
        for g in genome_files:
            print(genome_names[g], file=fp)

    np.set_printoptions(precision=3, suppress=True)
    print(mat)

if __name__ == '__main__':
    sys.exit(main())
