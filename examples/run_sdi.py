#!/usr/bin/env python

# import functions
import os, glob
import sdi


# autocorrelation to construct P reflectivity,
# at this step, the data would be saved in the SDI/ac path.
if True:
    filelst = glob.glob("SDI/good/IU.MBWA.00/*BHZ*SAC")
    sdi.autocorrelate(filelst)

# moveout correction
# at this step, the data would be saved in the SDI/mo0 path.
if True:
    filelst = glob.glob("SDI/ac/IU.MBWA.00/*BHZ*SAC")
    sdi.moveout(filelst, model="ak135", type="P", evdp_unit="m",
                twoside=True, slowness_ref=0.0)

# stacking over events
if True:
    pathlst = glob.glob("SDI/mo0/*")
    sdi.stack(pathlst=pathlst, extname="*.BHZ.M.SAC",
              taper_length=5.0, type="hann",
              filt=True, freqmin=0.5, freqmax=4, savepath="SDI/stk")

# AGC and time-to-depth conversion
if True:
    files = glob.glob("SDI/stk/*.sac")
    for file in files:
        print file
        sdi.data_for_plot(file, type="han", tlen=5.0, savepath="SDI/plot")

    # AGC
    files = glob.glob("SDI/plot/*.sac")
    for file in files:
        sdi.agc(file)

    # time-to-depth conversion using the ak135 model
    # the data in depth domain are named with the endings of ".d" in the path "SDI/plot/*.d"
    files = glob.glob("SDI/plot/*")
    for file in files:
        fn1, fn2 = os.path.splitext(file)
        if fn2 == ".d":
            continue
        sdi.time2depth_sac(file=file, dz=0.1, zrange=[0, 400])