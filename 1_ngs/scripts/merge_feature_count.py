#!/usr/bin/env python
import sys
import pandas as pd


def main():
    paths = sys.argv[1].split(",")
    samples = sys.argv[2].split(",")
    outfile = sys.argv[3]

    array = []
    for sample, path in zip(samples, paths):
        d = pd.read_csv(path, sep="\t", header=1, index_col=0)
        s = d[d.columns[-1]]
        s.name = sample
        array.append(s)
    dat = pd.concat(array, axis=1)
    dat = dat[sorted(dat.columns)]
    dat.to_csv(outfile, sep="\t")


if __name__ == "__main__":
    main()
