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
from sourmash.picklist import SignaturePicklist


def main():
    p = argparse.ArgumentParser()
    p.add_argument('sqldb')
    p.add_argument('idents', nargs='+')
    #p.add_argument('-o', '--output', required=True)
    args = p.parse_args()

    print(f"loading '{args.sqldb}'")
    db = LCA_SqliteDatabase(args.sqldb)
    ksize = db.ksize
    print(f"...done! got {len(db)} sketches at ksize={ksize}")

    # pick out just the requested sequences -
    picklist = SignaturePicklist('ident')
    picklist.init(values=args.idents)

    sub_db = db.select(picklist=picklist)

    #fp = open(args.output, 'w', newline="")
    #csv_w = csv.writer(fp)

    for n, sig in enumerate(sub_db.signatures()):
        if n % 100 == 0:
            print('...', n)

        mh = sig.minhash
        hashes = list(mh.hashes) # should be ordered, I think?
        max_size = len(hashes)

        hashval = hashes.pop()
        matches = set(db.hashval_to_idx[hashval])

        while hashes and len(matches) > 1:
            matches.intersection_update(db.hashval_to_idx[hashval])
            hashval = hashes.pop()

        print(hashes, len(matches))


if __name__ == '__main__':
    sys.exit(main())
