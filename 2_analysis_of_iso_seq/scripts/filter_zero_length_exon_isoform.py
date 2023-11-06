#!/usr/bin/env python
import sys

for line in sys.stdin:
    row = line.strip("\n").split("\t")
    flag = True
    for block_size in row[10].split(","):
        if int(block_size) <= 0:
            flag = False
            break
    if flag:
        sys.stdout.write(line)