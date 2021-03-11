#!/usr/bin/env python
import sys
import pandas as pd


def main():
    path1, path2 = sys.argv[1:]

    records = []
    with open(path1, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            columns = line.strip("\n").split("\t")
            columns[3] = int(columns[3]) - 1
            columns[4] = int(columns[4])
            attributes = dict()
            for item in columns[8].split(";"):
                item = item.strip()
                if item != "":
                    i = item.find("=")
                    assert i != -1
                    key = item[:i]
                    value = item[i + 1:]
                    value = value.strip("\"")
                    if key == "product":
                        value = value.replace("%2C", ",")
                    attributes[key] = value
            columns[-1] = attributes
            records.append(columns)

    header1 = ["Chrom", "Source", "Feature",
               "Start", "End", "Score", "Strand", "Frame"]
    header2 = []
    for record in records:
        for key in record[-1].keys():
            header2.append(key)
    header2 = list(sorted(set(header2)))
    header = header1 + header2
    with open(path2, "w+") as fw:
        fw.write("\t".join(header) + "\n")
        for record in records:
            values1 = record[:8]
            values2 = [record[-1].get(key, "") for key in header2]
            values = values1 + values2
            fw.write("\t".join(map(str, values)) + "\n")


if __name__ == "__main__":
    main()
