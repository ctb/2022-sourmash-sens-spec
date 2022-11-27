#! /usr/bin/env python
import sourmash
import sys
import argparse
import screed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('genomes', nargs='+')
    p.add_argument('-l', '--length', default=10000, type=int)
    args = p.parse_args()

    for filename in args.genomes:
        name = "_".join(filename.split('_')[:2])
        total_bp = 0
        filter_bp = 0
        for record in screed.open(filename):
            total_bp += len(record.sequence)
            if len(record.sequence) >= args.length:
                filter_bp += len(record.sequence)

        percent = f"{filter_bp / total_bp * 100:.1f}%"

        print(f"{name}, bp {total_bp:g}, filtered {filter_bp:g} {percent} (l={args.length:g})")


if __name__ == '__main__':
    sys.exit(main())
