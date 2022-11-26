#! /usr/bin/env python
"""
How many species have no perfectly informative hashvals?

Should be equivalent to the number of species where there are no genomes
with perfectly informative hashvals.
"""
import sys
import sourmash
from sourmash.lca import lca_utils
from sourmash import sourmash_args
from collections import defaultdict, Counter
import math
import pandas as pd
import argparse
import pickle


def identify_observed(db, hashvals, tax_observed):
    x = []

    for n, hashval in enumerate(hashvals):
        if n % 100000 == 0:
            print(f'... {n}')
            #if n: break

        lineages = db.get_lineage_assignments(hashval)
        tree = lca_utils.build_tree(lineages)
        lca, _ = lca_utils.find_lca(tree)
        if lca and lca[-1].rank == 'species':
            tax_observed[lca] = True

        #while lca:
        #    tax_observed[lca] = True
        #    lca = lca[:-1]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('db')
    args = p.parse_args()

    print(f"Loading '{args.db}'...")
    db = sourmash.load_file_as_index(args.db)
    print(f"...got {len(db)} sketches.")

    # build dictionary of all taxonomies at all levels
    tax_observed = {}
    for tax in db.lid_to_lineage.values():
        tax_observed[tax] = False
    print(f"Found {len(tax_observed)} total species.")

    identify_observed(db, db.hashvals, tax_observed)

    count = 0
    for k, v in tax_observed.items():
        if not v:
            count += 1

    print(f"missing {count} species.")

    if 0:
        n_rank_found = defaultdict(int)
        n_rank_total = defaultdict(int)
        for tax, is_found in tax_observed.items():
            while tax:
                rank = tax[-1].rank
                if is_found:
                    n_rank_found[rank] += 1
                n_rank_total[rank] += 1

                tax = tax[:-1]

        for k, total_count in n_rank_total.items():
            found_count = n_rank_found.get(k, 0)
            diff = total_count - found_count
            print(f"rank {k}, found: {found_count}, total: {total_count}, diff: {diff}")

    
if __name__ == '__main__':
    sys.exit(main())
