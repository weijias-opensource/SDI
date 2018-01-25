#!/usr/bin/env python

"""
Codes for calculation the instantaneous frequency of a real signal using hilbert transform
"""


import numpy as np
from scipy.signal import hilbert, savgol_filter

def envelope(xa):
    """
    The amplitude envelope: magnitude of the analytical signal

    Parameters
    ----------
    xa : array_like
         analytical signal

    Returns
    -------
    env : envelope of the analytical signal
    """
    return np.abs(xa)

def instantaneous_phase(xa):
    """
    The instantaneous phase

    Parameters
    xa : analytical signal

    Returns
    phase: instantaneous phase
    """
    return np.unwrap(np.angle(xa))

def instantaneous_frequency(xa, fs):
    """
    The instantaneous frequency

    :param xa:
    :param fs: frequency of sampling
    :return frequency:
    """
    phase = instantaneous_phase(xa)
    return np.diff(phase) / (2.0*np.pi)*fs

def instantaneous_frequency_smooth(xa, fs):
    """
    The instantaneous frequency smoothed locally

    :param xa:
    :param fs: frequency of sampling
    :return frequency:
    """
    phase = instantaneous_phase(xa)
    phase_dev = savgol_filter(phase, window_length=7, polyorder=3, deriv=1)

    # freq = np.diff(phase2) / (2.0*np.pi)*fs
    # a = np.zeros((1,), dtype=freq.dtype)

    # print freq.shape, a.shape
    # freq = np.concatenate((freq, a))

    # print len(freq), len(phase), len(xa)

    freq = phase_dev / (2.0*np.pi)*fs * np.conj(xa) / (xa * np.conj(xa))
    return freq



if __name__=="__main__":
    from obspy import read
    import matplotlib.pyplot as plt


    st = read("2016.092.02.39.05.0000.ZJ.ZHS.00.BHZ.M.SAC")
    tr = st[0]
    fs = 1.0 / tr.stats.delta

    x = tr.data

    print fs

    print len(x)

    xa = hilbert(x)
    freq = instantaneous_frequency_smooth(xa, fs)
    # ph, ph2 = instantaneous_frequency_smooth(xa, fs)

    tr.data = freq

    tr.write("test_smooth.sac", format="SAC")

    # fig = plt.figure()
    # ax0 = fig.add_subplot(211)
    # ax0.plot(x)
    # ax1 = fig.add_subplot(212)
    # ax1.plot(freq)
    # # ax1.plot(ph)
    # # ax1.plot(ph2,"r")
    #
    # plt.show()
