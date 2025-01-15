#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

import pod5 as p5
import sys

# Provide pod5 file and read ID
pod5_file = str(sys.argv[1])
selected_read_id = str(sys.argv[2])

with p5.Reader(pod5_file) as reader:

    # Read the selected read from the pod5 file
    # next() is required here as Reader.reads() returns a Generator
    read = next(reader.reads(selection=[selected_read_id]))

    # Get the signal data and sample rate
    sample_rate = read.run_info.sample_rate
    signal = read.signal

    # Compute the time steps over the sampling period
    time = np.arange(len(signal)) / sample_rate

    plt.figure().set_figwidth(15)
   #plt.xlim(0, 7)
    plt.ylim(200, 1200)

    # Plot using matplotlib
    plt.plot(time, signal)

    # show the plot
    plt.savefig(selected_read_id + '.jpg')
