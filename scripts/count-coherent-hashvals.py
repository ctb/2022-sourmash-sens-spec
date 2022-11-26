#! /usr/bin/env python
import sys
import sourmash
from sourmash.lca import lca_utils
from sourmash import sourmash_args
from collections import defaultdict, Counter
import math
import pandas as pd
import argparse
import pickle


def calc_entropy(lineages, *, rank='species'):
    cnt = Counter()
    for lin in lineages:
        lin = lca_utils.pop_to_rank(lin, rank)
        cnt[lin] += 1

    assert cnt, lineages

    total = sum(cnt.values())
    H = 0
    for v in cnt.values():
        p = v / total
        H -= p * math.log(p, 2)

    return H


def entropy_of_lineages(db, hashvals, *, rank='species'):
    x = []

    for n, hashval in enumerate(hashvals):
        if n % 100000 == 0:
            print(f'... {n}')
            #if n: break

        lineages = db.get_lineage_assignments(hashval)
        H = calc_entropy(lineages, rank=rank)
        x.append(dict(hashval=hashval, H=H))

    return pd.DataFrame(x)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('db')
    p.add_argument('-o', '--output-signature', required=True,
                   help='signature containing above-family hashvals')
    args = p.parse_args()

    print(f"Loading '{args.db}'...")
    db = sourmash.load_file_as_index(args.db)
    print(f"...got {len(db)} sketches.")

    species_df = entropy_of_lineages(db, db.hashvals, rank='species')

    num_species_h0 = len(species_df[species_df.H == 0.0])
    species_nonzero_df = species_df[species_df.H > 0.0]
    total = len(species_df)

    print(f"{num_species_h0} of {total} hashvals ({num_species_h0 / total * 100:.1f}%) are perfectly informative at species level!")

    print(f"examining remaining {len(species_nonzero_df)} at genus level.")
    genus_df = entropy_of_lineages(db, species_nonzero_df.hashval, rank='genus')

    num_genus_h0 = len(genus_df[genus_df.H == 0.0])
    genus_nonzero_df = genus_df[genus_df.H > 0.0]

    print(f"{num_genus_h0} of {total} hashvals ({num_genus_h0 / total * 100:.1f}%) are perfectly informative at genus level!")

    print(f"examining remaining {len(genus_nonzero_df)} at family level.")
    family_df = entropy_of_lineages(db, genus_nonzero_df.hashval, rank='family')

    num_family_h0 = len(family_df[family_df.H == 0.0])
    family_nonzero_df = family_df[family_df.H > 0.0]

    print(f"{num_family_h0} of {total} hashvals ({num_family_h0 / total * 100:.1f}%) are perfectly informative at family level!")

    # to save, make an empty copy of the MinHash used in this database -
    first_ss = next(iter(db.signatures()))
    first_mh = first_ss.minhash

    mh = first_mh.copy_and_clear()
    
    hashvals = set(family_nonzero_df.hashval)
    mh.add_many(hashvals)

    output_ss = sourmash.SourmashSignature(mh)
    with sourmash_args.SaveSignaturesToLocation(args.output_signature) as save_sig:
        save_sig.add(output_ss)
    print(f"Saved {len(hashvals)} hashvals to '{args.output_signature}'")

    
if __name__ == '__main__':
    sys.exit(main())
