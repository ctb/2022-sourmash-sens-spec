#! /usr/bin/env python
"""
defining unicity here as "the minimum number of hashes necessary to
uniquely identify this genome in GTDB."
"""
import sys
import csv
import argparse
from collections import defaultdict
import random

import sourmash
from sourmash.index.sqlite_index import LCA_SqliteDatabase


def main():
    p = argparse.ArgumentParser()
    p.add_argument('sqldb')
    p.add_argument('-o', '--output', required=True)
    args = p.parse_args()

    print(f"loading '{args.sqldb}'")
    db = LCA_SqliteDatabase(args.sqldb)
    ksize = db.ksize
    print(f"...done! got {len(db)} sketches at ksize={ksize}")

    fp = open(args.output, 'w', newline="")
    csv_w = csv.writer(fp)

    for n, sig in enumerate(db.signatures()):
        if n % 100 == 0:
            print('...', n)

        mh = sig.minhash
        hashes = list(mh.hashes) # should be ordered, I think

        # probably unnecessary?
        #random.shuffle(hashes)

        hashval = hashes.pop()
        matches = set(db.hashval_to_idx[hashval])

        matches_at = [len(matches)]

        while hashes and len(matches) > 1:
            matches.intersection_update(db.hashval_to_idx[hashval])
            matches_at.append(len(matches))
            hashval = hashes.pop()

        name = sig.name.split(' ')[0]
        csv_w.writerow([name, sig.md5sum()] + matches_at)


if __name__ == '__main__':
    sys.exit(main())
