#!/usr/bin/env python

import sys
import pandas as pd


def main():
    path1, path2, path3 = sys.argv[1:]

    dat1 = pd.read_csv(path1, sep="\t", index_col=0, header=0)
    dat2 = pd.read_csv(path2, sep="\t", index_col=0, header=0)
    columns = list(dat2.columns)
    columns.remove("Length")
    dat2 = dat2[columns]
    dat = dat1.merge(dat2, left_index=True, right_index=True)
    dat.to_csv(path3, sep="\t")


if __name__ == "__main__":
    main()