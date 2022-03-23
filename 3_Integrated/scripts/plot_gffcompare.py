#!/usr/bin/env python
import sys
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams["font.family"] = "arial"

infile, outfile = sys.argv[1:]
dat = pd.read_csv(infile, sep="\t", header=None)
counter = Counter(dat[3])
items = list(sorted(counter.items(), key=lambda item: item[1], reverse=True))
xticks = [item[0] for item in items]
counts = [item[1] for item in items]
xs = np.arange(len(xticks))
ys = np.array(counts) * 100 / sum(counts)
plt.figure(figsize=(4, 2.5))
plt.bar(xs, ys, color="black")
plt.xlim(xs[0] - 0.5, xs[-1] + 0.5)
plt.ylim(-1, 40)
plt.xlabel("Category")
plt.ylabel("Percentage (%)")
plt.xticks(xs, xticks)
plt.tight_layout()
plt.savefig(outfile, dpi=300)
plt.close()

