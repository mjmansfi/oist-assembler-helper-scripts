#!/usr/bin/env python3

# remember to mount the working directory so python within the container can access it
import sys
from scipy.signal import argrelmin
from scipy.signal import find_peaks
import logging
import numpy as np

# Log configuration
logging.basicConfig(filename='purge_haplotigs_minima.log', level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.debug("Starting script...")

hist_dict = {}

logging.debug("Reading %s..." % sys.argv[1])
with open(sys.argv[1]) as i:
    for line in i:
        line = line.strip()
        # Keys are sequencing coverage values, values are counts of that coverage value
        hist_dict[line.split(',')[0]] = int(line.split(',')[1])

"""
In an unphased genome assembly, the haplotigs will exhibit around
half of the coverage of the diploid co-assembled contigs. This can
be exploited to separate haplotigs: see the tutorial in

https://bitbucket.org/mroachawri/purge_haplotigs/wiki/Tutorial

Under the assumption that the input assembly exhibits the expected
pattern of coverage, the low, midpoint, and high cutoffs required by
phase_haplotigs can be estimated as the local minima of the coverage
graph when the order of the coverage vs. count equals 3.

This assumption cannot be true in all situations. Check the graphical
output of purge_haplotigs readhist.
"""
coverage_data = np.array(list(hist_dict.values()))
num_peaks = len(find_peaks(coverage_data, width=10)[0])
num_peaks_message = "%s peaks were estimated from the histogram generated from your .bam file." % num_peaks
logging.info(num_peaks_message)

if num_peaks < 2:
    num_peaks_warning = "WARNING:\t%s coverage peaks were estimated in the histogram generated from your .bam file,\n\twhen 2 or more are expected. The low, midpoint, and high cutoff points estimated\n\tby this script may not be sensible. Please inspect the pattern in the purge_haplotigs .png output.\n\tFor more information, see:\n\thttps://bitbucket.org/mroachawri/purge_haplotigs/wiki/Tutorial\n" % num_peaks
    logging.warning(num_peaks_warning)

local_minima = argrelmin(coverage_data, order = 3)
local_minima_x = [list(hist_dict.keys())[X] for X in local_minima[0]]
local_minima_y = [list(hist_dict.values())[X] for X in local_minima[0]]

with open("critical_values.csv", "w") as cv_out:
    cv_out.write("Critical point, X, Y\nLow, %s, %s\nMidpoint, %s, %s\nHigh, %s, %s\n" % (local_minima_x[0], local_minima_y[0], local_minima_x[1], local_minima_y[1], local_minima_x[2], local_minima_y[2]))

with open("low_mid_high.csv", "w") as lmh_out:
    lmh_out.write("cutoff_low,cutoff_mid,cutoff_high\n%s,%s,%s\n" % (local_minima_x[0], local_minima_x[1], local_minima_x[2]))
