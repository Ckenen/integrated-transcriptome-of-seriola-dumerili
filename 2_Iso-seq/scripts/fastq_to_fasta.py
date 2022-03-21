#!/usr/bin/env python

import sys


def main():
    path1 = sys.argv[1]
    path2 = sys.argv[2]

    with open(path1) as f, open(path2, "w+") as fw:
        for i, line in enumerate(f):
            j = i % 4
            if j == 0:
                line = ">" + line[1:]
                fw.write(line)
            elif j == 1:
                fw.write(line)

if __name__ == "__main__":
    main()
