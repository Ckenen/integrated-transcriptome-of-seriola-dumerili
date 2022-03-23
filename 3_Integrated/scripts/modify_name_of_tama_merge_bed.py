#!/usr/bin/env python

import sys


def main():
    path1, path2 = sys.argv[1:]

    with open(path1) as f, open(path2, "w+") as fw:
        for line in f:
            cols = line.strip("\n").split("\t")
            cols[3] = ";".join(cols[3].split(";")[:2])
            cols[6] = cols[1]
            cols[7] = cols[1]
            fw.write("\t".join(cols) + "\n")
            

if __name__ == "__main__":
    main()