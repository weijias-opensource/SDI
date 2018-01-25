#!/usr/bin/env python

# import functions
import os, glob
from obspy import read
import matplotlib.pyplot as plt
import numpy as np

files = glob.glob("SDI/plot/*.d")
print "number of files", len(files)

files.sort()

tr1 = read(files[0])[0]
tr2 = read(files[1])[0]

depth = np.arange(len(tr1.data)) * tr1.stats.delta

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(4,6))
ax1.plot(tr1.data, depth, linewidth=0.5)
ax2.plot(tr2.data, depth, linewidth=0.5)

ax1.set_ylabel("Depth [km]")
ax1.set_ylim([300, 0])
ax2.set_ylim([300, 0])
ax1.set_title("AGC")
ax1.set_title("No AGC")
plt.tight_layout()

plt.show()