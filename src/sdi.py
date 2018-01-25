#!/usr/bin/env python


import glob
import os
import math

import hilbert

import numpy as np

from obspy import read
from obspy.taup import TauPyModel
from obspy.core.trace import Trace
from obspy.geodetics import locations2degrees

from scipy import signal
from scipy.interpolate import interp1d
# from obspy.signal.util import nextpow2
from scipy.fftpack import fft, ifft

# def walk():
#     pass

def preprocess(filelst, freqmin=0.05, freqmax=5.0, evdp_unit="m", sample_rate=40):
    """
    Perform preprocessing of data, including dedtrend, taper, filter and calculate
    theoretical travel time for P and S.

    The theoretical P and S traveltime would be saved in header t1 and t2.

    The processed SAC file would write into the SDI/data/NET.STA/

    :param filelst:
    :param freqmin: min fre, default 0.05
    :param freqmax: max fre, default 5.0
    :param evdp_unit: "km" or "m", default km
    :param sample_rate: interpolate the traces in the the given sampling rate, default 40
    :return: none
    """


    # try:
    #     os.makedirs("data")

    model = TauPyModel(model="ak135")

    if evdp_unit=="m":
        evdp = 1000.0
    elif evdp_unit=="km":
        evdp = 1.0
    else:
        evdp = 1.0

    filelst.sort()


    for file in filelst:

        print "Pre-processing:", file

        try:
            st = read(file)
            tr = st[0]
        except:
            continue

        # pre-processing data
        tr1 = tr.copy()
        tr1.detrend(type="linear")
        tr1.detrend(type="demean")
        tr1.taper(type="cosine", max_percentage=0.05, side="both")
        tr1.filter(type="bandpass", freqmin=freqmin, freqmax=freqmax, zerophase=True)

        if tr.stats.sampling_rate!=sample_rate:
            # tr1.interpolate(sampling_rate=sample_rate)
            tr1.resample(sampling_rate=sample_rate)

        # calculating the traveltime of P, S
        arrivals = model.get_travel_times(source_depth_in_km=tr1.stats.sac.evdp/evdp,
                                          distance_in_degree=tr1.stats.sac.gcarc,
                                          phase_list=["P", "S"])

        for arr in arrivals:
            if arr.name == "P":
                tr1.stats.sac.t1 = arr.time + tr1.stats.sac.o
                tr1.stats.sac.kt1 = "P"
                tr1.stats.sac.user1 = arr.ray_param*math.pi/180.0/111.195
                tr1.stats.sac.kuser1 = "P"
            elif arr.name =="S":
                tr1.stats.sac.t2 = arr.time + tr1.stats.sac.o
                tr1.stats.sac.kt2 = "S"
                tr1.stats.sac.user2 = arr.ray_param*math.pi/180.0/111.195
                tr1.stats.sac.kuser2 = "S"

        # event origin time
        evt_origin = tr1.stats.starttime + tr1.stats.sac.o - tr1.stats.sac.b
        o = str(evt_origin)
        for s in ["-", "T", ":"]:
            o = o.replace(s, "_")

        # fn = "Event_"+o[:19]+"."+tr1.stats.channel+".sac"

        temp = evt_origin.datetime
        fn = ".".join([temp.strftime("%Y.%j.%H.%M.%S"), "0000", tr1.id, "M", "SAC"])

        # print o

        # save sac files
        folder = "SDI/data/" + ".".join([tr1.stats.network, tr1.stats.station])

        if tr1.stats.location!="":
            folder = folder+"."+tr1.stats.location

        # print folder


        try:
            os.makedirs(folder)
        except:
            pass

        # print folder

        # print tr1.stats.channel
        fn = folder+"/"+fn
        tr1.write(fn, format="SAC")

        # return

    pass

def snr(filelst, threshold=2.0, sigmin=-10.0, sigmax=10.0, noimin=-100.0, noimax=-10.0, id="P", cutmin=-30, cutmax=150):
    """
    S/N ratio calculation and choose traces with SNR>threshold.

    The signal is defined in the time window sigmin seconds and sigmax seconds relative to the theoreticla travel time.
    The noise is defined in the time window noimin s and noimax s relative to the theoretical traveltime.
    The retrieved waveform is in the time windows between cutmin and cutmax relative to the theoretical traveltime.

    The chosen files would be saved in the path SDI/good/NET.STA


    :param filelst:
    :param threshold: The traces with SNR greater than threshold would be chosen, default 2.0
    :param sigmin: min signal time window, default -10
    :param sigmax: max signal time window
    :param noimin: min noise time window
    :param noimax: max noise time window
    :param id: "P" or "S", default P
    :param cutmin: min time window of extracted waveform
    :param cutmax: max time window of extracted waveform
    :return:
    """

    for file in filelst:

        tr = read(file)[0]

        if id=="P":
            try:
                tref = tr.stats.sac.t1
            except:
                continue
        elif id == "S":
            try:
                tref = tr.stats.sac.t2
            except:
                continue
        else:
            print "ERROR: id must be P or S."
            os._exit(0)

        if tref=="-12345.0":
            continue

        # print tref, tr.stats.starttime, tr.stats.sac.b, tr.stats.sac.o

        tref_abs = tr.stats.starttime + tref - tr.stats.sac.b

        # print tref_abs

        # calculating SNR
        sigstart = tref_abs + sigmin
        sigend   = tref_abs + sigmax

        noistart = tref_abs + noimin
        noiend   = tref_abs + noimax

        tr1 = tr.slice(starttime=sigstart, endtime=sigend)
        tr2 = tr.slice(starttime=noistart, endtime=noiend)

        sig = tr1.data
        noi = tr2.data

        # print sig.shape, noi.shape

        sig = np.sqrt(np.sum(sig*sig))/len(sig)
        noi = np.sqrt(np.sum(noi*noi))/len(noi)

        # print "SNR=", sig/noi, tr.stats.channel

        if noi<=0:
            continue

        try:
            snr = sig/noi
        except:
            continue

        if snr<threshold:
            continue

        tstart = tref_abs + cutmin
        tend   = tref_abs + cutmax
        tr3 = tr.slice(starttime=tstart, endtime=tend)


        path, fn = filen2(file)
        fullfn = "SDI/good/"+path+"/"+fn

        # print fullfn
        try:
            os.makedirs("SDI/good/"+path)
        except:
            pass

        tr3.write(fullfn, format="SAC")

        print "SNR: ", snr, file



    pass

def filen2(sacfile):
    """
    new version of filen. Return path and filename

    :param sacfile:
    :return:
    """


    tr1 = read(sacfile)[0]

    # event origin time
    evt_origin = tr1.stats.starttime + tr1.stats.sac.o - tr1.stats.sac.b
    o = str(evt_origin)
    for s in ["-", "T", ":"]:
        o = o.replace(s, "_")

    # filename = "Event_"+o[:19]+"."+tr1.stats.channel+".sac"
    temp = evt_origin.datetime
    filename = ".".join([temp.strftime("%Y.%j.%H.%M.%S"), "0000", tr1.id, "M", "SAC"])

    # print fn

    # save sac files
    folder = ".".join([tr1.stats.network, tr1.stats.station])
    if tr1.stats.location!="":
        folder = folder+"."+tr1.stats.location

    path = folder
    return path, filename

    pass

def filen(sacfile):
    """
    make a filename from sac file using station
    :param sacfile: input filename of SAC file
    :return: path, network.station.loc (if loc not available, then network.station)
             filename, Event_2016_10_01_10_01_23_BHE.sac

    """

    tr1 = read(sacfile)[0]

    # event origin time
    evt_origin = tr1.stats.starttime + tr1.stats.sac.o - tr1.stats.sac.b
    o = str(evt_origin)
    for s in ["-", "T", ":"]:
        o = o.replace(s, "_")

    filename = "Event_"+o[:19]+"."+tr1.stats.channel+".sac"

    # print fn

    # save sac files
    folder = ".".join([tr1.stats.network, tr1.stats.station])
    if tr1.stats.location!="":
        folder = folder+"."+tr1.stats.location

    path = folder
    return path, filename


    # try:
    #     os.makedirs(folder)
    # except:
    #     pass

def autocorrelate(filelst):
    """
    Autocorrelation in frequency domain.


    The implementation of AC is from SCIPY. The autocorrelograms is normalized.

    :param filelst:
    :return:
    """

    for file in  filelst:

        print "AC: ", file

        tr = read(file)[0]
        x = tr.data

        stel = tr.stats.sac.stel
        stla = tr.stats.sac.stla
        stlo = tr.stats.sac.stlo

        try:
            # auto-correlation from SCIPY implemented in F domain, which is much faster than those in time domain
            ac = signal.fftconvolve(x, x[::-1], mode="full")
            ac = ac/max(abs(ac))
        except:
            continue

        tr1 = Trace(data=ac)
        tr1.stats.network = tr.stats.network
        tr1.stats.station =tr.stats.station
        tr1.stats.location = tr.stats.location
        tr1.stats.delta = tr.stats.delta
        tr1.stats.sampling_rate = tr.stats.sampling_rate
        tr1.stats.channel = tr.stats.channel

        tr1.stats._format = "SAC"
        tr1.stats.sac = {u'stla': stla,  u'stlo': stlo, u'stel': stel}
        tr1.stats.sac.user1 = tr.stats.sac.user1
        tr1.stats.sac.user2 = tr.stats.sac.user2
        tr1.stats.sac.kuser1 = tr.stats.sac.kuser1
        tr1.stats.sac.kuser2 = tr.stats.sac.kuser2

        # tr1.stats.sac.stla = tr.stats.sac.stla
        # tr1.stats.sac.stlo = tr.stats.sac.stlo
        # tr1.stats.sac.stel = tr.stats.sac.stel


        path, fn = filen2(file)
        fullfn = "SDI/ac/"+path+"/"+fn

        # print fullfn
        try:
            os.makedirs("SDI/ac/"+path)
        except:
            pass

        tr1.write(fullfn, format="SAC")


    pass


def autocorrelate2(filelst):
    """
    Autocorrelation in frequency domain.


    The implementation of AC is from SCIPY. The autocorrelograms is normalized.

    :param filelst:
    :return:
    """

    for file in  filelst:

        print "AC: ", file

        tr = read(file)[0]
        x = tr.data

        stel = tr.stats.sac.stel
        stla = tr.stats.sac.stla
        stlo = tr.stats.sac.stlo

        try:
            # auto-correlation from SCIPY implemented in F domain, which is much faster than those in time domain
            ac = signal.fftconvolve(x, x[::-1], mode="full")
            ac = ac/max(abs(ac))
        except:
            continue

        tr1 = Trace(data=ac)
        tr1.stats.network = tr.stats.network
        tr1.stats.station =tr.stats.station
        tr1.stats.location = tr.stats.location
        tr1.stats.delta = tr.stats.delta
        tr1.stats.sampling_rate = tr.stats.sampling_rate
        tr1.stats.channel = tr.stats.channel

        tr1.stats._format = "SAC"
        tr1.stats.sac = {u'stla': stla,  u'stlo': stlo, u'stel': stel}
        tr1.stats.sac.user1 = tr.stats.sac.user1
        tr1.stats.sac.user2 = tr.stats.sac.user2
        tr1.stats.sac.kuser1 = tr.stats.sac.kuser1
        tr1.stats.sac.kuser2 = tr.stats.sac.kuser2

        # tr1.stats.sac.stla = tr.stats.sac.stla
        # tr1.stats.sac.stlo = tr.stats.sac.stlo
        # tr1.stats.sac.stel = tr.stats.sac.stel


        path, fn = filen2(file)
        fullfn = "SDI/ac/"+path+"/"+fn

        # print fullfn
        try:
            os.makedirs("SDI/ac/"+path)
        except:
            pass

        tr1.write(fullfn, format="SAC")


    pass

def stack(pathlst, extname="*.BHZ.sac", taper_length=5.0, type="hann", filt=True, freqmin=0.5, freqmax=4, savepath="SDI/stk"):
    """
    Stacking autocorrelograms over each stations.

    :param pathlst:
    :param extname: used for glob files to  be stacked
    :param taper_length: seconds, default 5.0
    :param type: "hann" or "zero"
    :param filt: True or False, default True
    :param freqmin: min freq
    :param freqmax: max freq
    :param savepath: where to save stacked data
    :return:
    """

    try:
        # os.makedirs("SDI/stk")
        os.makedirs(savepath)
    except:
        pass

    pathlst.sort()



    for path in pathlst:
        filelst = glob.glob(path+"/"+extname)

        file = filelst[0]
        tr = read(file)[0]

        # remove impulse at t=0
        delta = tr.stats.delta
        npts = tr.stats.npts
        ntaper = int(round(taper_length/delta))*2

        if type=="hann":
            hann = np.hanning(ntaper)
            hann = hann*hann
        elif type=="zero":
            hann = np.zeros(ntaper)

        a1 = hann[:ntaper/2]
        a2 = hann[ntaper/2:]
        hann = np.ones_like(tr.data)
        a1 = a1[::-1]
        a2 = a1[::-1]
        i1 = npts/2 - len(a1)
        i2 = npts/2
        hann[i1:i2] = a1
        i1 = npts/2
        i2 = npts/2 + len(a2)
        hann[i1:i2] = a2

        summ = np.zeros_like(tr.data)
        try:
            count = 0
            for file in filelst:
                tr = read(file)[0]
                if tr.stats.npts!=npts:
                    continue
                summ += tr.data
                count += 1
        except:
            continue

        tr.data=summ
        if filt:
            tr.filter(type="bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
        tr.normalize()
        tr.data = tr.data * hann


        fn = savepath+"/"+tr.stats.network+"."+tr.stats.station
        if tr.stats.location!="":
            fn += "."+tr.stats.location

        fn += "."+tr.stats.channel+".sac"

        print "STK: ", fn, count
        tr.write(filename=fn, format="SAC")


    pass

# reference Kennett et al., 2015 GRL, 065345.
def spatial_stack(filelst, sigma=0.5, weight_threshold=1.2e-4, filt=True, freqmin=0.5, freqmax=4, savepath="SDI/spatial_stack"):

    # savepath = "SDI/spatial_stack"
    try:
        os.makedirs(savepath)
    except:
        pass

    filelst.sort()
    for file1 in filelst:
        # print "spatial stacking:", file1

        basename = os.path.basename(file1)

        tr1 = read(file1)[0]
        summ = np.zeros_like(tr1.data)

        lat1 = tr1.stats.sac.stla
        lon1 = tr1.stats.sac.stlo

        sum_weight = 0.0

        for file2 in filelst:
            tr2 = read(file2)[0]

            if tr1.stats.channel!=tr2.stats.channel:
                continue

            lat2 = tr2.stats.sac.stla
            lon2 = tr2.stats.sac.stlo

            dist = locations2degrees(lat1=lat1, long1=lon1, lat2=lat2, long2=lon2)

            weight = math.exp(-dist*dist/sigma/sigma)
            if weight<weight_threshold:
                continue

            if filt:
                tr2.filter(type="bandpass", freqmin=freqmin, freqmax=freqmax)

            tr2.normalize()

            try:
                sum_weight += weight
                summ += tr2.data*weight
            except:
                continue

            # print file2, dist, weight

        # print sum_weight

        # if sum_weight<1.5:
        #     print sum_weight
        #     continue

        tr1.data = summ/sum_weight

        tr1.write(filename=savepath+"/"+basename, format="SAC")
        print savepath+"/"+basename


        # return


    pass


##############################################################################
# The following are codes of move out
##############################################################################
def time2depth_test(model="ak135", twoway=True):
    """
    convert from time to depth
    twoway, two-way time
    """

    # model in time-domain and depth domain
    depth, vp, vs, tvp, tvs = d2t(model=model, twoway=twoway)

    t = np.arange(0.0, 120.025, 1.0)

    f = interp1d(tvp, depth)
    d = f(t) # depth for P

    fs = interp1d(tvs, depth)
    ds = fs(t) # depth for S

    for i in range(len(t)):
        print "%f\t%f\t%f" % (t[i], d[i], ds[i])


def time2depth(time, type="P", model="ak135"):

    depth, vp, vs, tvp, tvs = d2t(model=model, twoway=True)

    f = interp1d(tvp, depth)
    d = f(time)

    return d

    pass

def d2t(model="ak135", twoway=True):

    with open(model + ".dat", "r") as fp:
        lst = fp.readlines()

    depth = []
    vp = []
    vs = []
    ro = []

    for line in lst:
        row = line.split()

        v1 = float(row[0])
        v2 = float(row[1])
        v3 = float(row[2])
        v4 = float(row[3])

        if v1 not in depth:
            depth.append(v1)
            vp.append(v2)
            vs.append(v3)
            ro.append(v4)

    depth = np.array(depth)
    vp = np.array(vp)
    vs = np.array(vs)
    ro = np.array(ro)

    tvp = np.zeros_like(depth)
    tvs = np.zeros_like(depth)

    for i in range(1, len(depth)):
        tvp[i] = (depth[i] - depth[i - 1]) / vp[i]
        tvs[i] = (depth[i] - depth[i - 1]) / vs[i]

    # tvp[0] = tvp[1]
    # tvs[0] = tvs[1]

    for i in range(1, len(depth)):
        tvp[i] += tvp[i - 1]
        tvs[i] += tvs[i - 1]

    if twoway:  # two-way travel time corresponding to the ak135 model
        tvp = tvp * 2.0
        tvs = tvs * 2.0

    # print tvp

    return depth, vp, vs, tvp, tvs

    pass


def model_in_time_domain(delta, npts, type="P", model="ak135", twoway=True):
    """
    convert velocity model in depth domain into time domain
    return: tnew same as the time series of reflections
            v the corresponding velocity of tnew, the mean velocity.
    """
    depth, vp, vs, tvp, tvs = d2t(model=model, twoway=True)  # vp vs interval velocity

    vpa = np.zeros_like(vp)
    vsa = np.zeros_like(vs)

    tvp = np.array(tvp)
    tvs = np.array(tvs)

    # tvp[0] = 0.01
    # tvs[0] = 0.01

    vpa = 2.0 * depth / tvp
    vsa = 2.0 * depth / tvs
    vpa[0] = vpa[1]
    vsa[0] = vsa[1]

    # tvp[0] = 0
    # tvs[0] = 0

    tnew = np.arange(npts) * delta

    if type == "P" or type == "p":
        f = interp1d(tvp, vpa)
        v = f(tnew)
    elif type == "S" or type == "s":
        f = interp1d(tvs, vsa)
        v = f(tnew)
    else:
        print "model_in_time_domain: type must be P or S."
        os._exit(0)

    return tnew, v


def moveout(filelst, model="ak135", type="P", evdp_unit="m", twoside=True, slowness_ref=6.0):
    """
    moveout according to the first order relationship
    profile output time profile or depth profile
    p_ref = 6.0 for P unit: s/deg
    """

    slowness_ref /= 111.195
    # print slowness_ref
    # return

    for file in filelst:

        print "Moveout: ", file

        tr = read(file)[0]

        npts = tr.stats.npts/2
        delta = tr.stats.delta
        data = tr.data[npts:]
        npts = len(data)

        if type == "P" or type== "p":
            slo = tr.stats.sac.user1
            tr.stats.sac.user1 = slowness_ref
        elif type == "S" or type== "s":
            slo = tr.stats.sac.user2
            tr.stats.sac.user2 = slowness_ref

        # print slo
        if slo < 0:
            print file, ": phase is not valid."
            continue

        # moveout
        # velocity model in time domain
        t, v = model_in_time_domain(delta, npts, type=type, model=model, twoway=True)


        # res = 1.0 - 0.5 * v * v * slo * slo  # B.L.N. Kennett gives the minus rather than plus. Talk it later.
        # res = np.sqrt(1.0 - v * v * slo * slo)
        # res = np.sqrt(1.0 - v * v * slo * slo) #- np.sqrt(1.0 - v * v * slowness_ref * slowness_ref)
        res = 1.0 - 0.5 * v * v *(slo*slo - slowness_ref*slowness_ref)
        tau0 = t / res

        d_interp = np.zeros_like(t)

        ft = interp1d(tau0, data)

        t_temp = np.arange(0, max(tau0), delta)
        d_temp = np.zeros_like(t_temp)

        d_temp = ft(t_temp)

        if len(t_temp) < len(t):
            d_interp[0:len(d_temp)] = d_temp
        else:
            d_interp = d_temp[0:len(d_interp)]


        if twoside:
            d_sys = d_interp[1:]
            d_sys = d_sys[::-1]
            d_interp = np.concatenate((d_sys, d_interp))


        tr.data = d_interp

        # save data
        # since the event information is erased in autocorrelation, so we use the basename here
        localpath = tr.stats.network+"."+tr.stats.station
        if tr.stats.location!="":
            localpath += "." + tr.stats.location
        fn = os.path.basename(file)

        if slowness_ref==0.0:
            fullfn = "SDI/mo0/"+localpath+"/"+fn
        else:
            fullfn = "SDI/mo/"+localpath+"/"+fn

        try:
            if slowness_ref==0.0:
                os.makedirs("SDI/mo0/"+localpath)
            else:
                os.makedirs("SDI/mo/"+localpath)
        except:
            pass

        tr.write(fullfn, format="SAC")



    pass


##############################################################################
# end moveout
##############################################################################

def data_for_plot(file, type="han", tlen=5.0, savepath="SDI/plot"):
    """
    Preparing data for plotting. Input data should be symmetric.

    Taper the first tlen second

    :param file:
    :param type: "han" or "zero"
    :param tlen: taper window length in seconds
    :param savepath: save path for data
    :return:
    """

    # read data from file
    tr = read(file)[0]

    npts = tr.stats.npts
    nhalf = npts/2
    tr.data = tr.data[nhalf:]
    tr.stats.npts = len(tr.data)

    times = 3.0
    if type == "han":
        hann = np.ones_like(tr.data)
        ntaper = int(round(tlen / tr.stats.delta))
        hann1 = np.hanning(ntaper*2)
        hann[:ntaper] = hann1[:ntaper]
        hann = np.power(hann, times)
        tr.data *= hann

    elif type == "zero":
        nt = int(tlen/tr.stats.delta)
        tr.data[0:nt] = 0.0

    basenm = os.path.basename(file)
    fullname = savepath + "/" + basenm

    try:
        os.makedirs(savepath)
    except:
        pass

    tr.write(filename=fullname, format="SAC")



def ac_hilbert(file, envelope=True):
    """

    Calculating Instantaneous Frequency

    :param file:
    :param envelope: If True, then output envelope of AC and normalized AC by envenlope
    :return:
    """

    tr = read(file)[0]
    delta = tr.stats.delta
    data = tr.data

    fs = 1 / delta
    xa = hilbert.hilbert(tr.data)
    freq = hilbert.instantaneous_frequency(xa, fs)
    tr.data = freq

    # for i in range(0, len(tr.data)):
    #     d = tr.data[i]
    #     if d > 0.0:
    #         tr.data[i] = 0.0

    fn, fext = os.path.splitext(file)

    if fext == ".norm":
        myext = ".ifn"
    else:
        myext = ".if"

    fullname = fn + myext
    print fullname

    tr.write(filename=fullname, format="SAC")

    if envelope:
        fullname = fn + ".env"
        tr.data = hilbert.envelope(xa)
        tr.write(filename=fullname, format="SAC")

        fullname = fn + ".norm"
        tr.data = data / tr.data
        tr.data = np.nan_to_num(tr.data)
        tr.write(filename=fullname, format="SAC")

    pass

def spectral_whiten(file, minfreq=0.5, maxfreq=4.0):
    pass

def agc(file):

    tr = read(file)[0]
    d = tr.data

    # scale = 2.0
     # rms = np.sqrt(np.sum(d*d/len(d))) * scale

    fn, fext = os.path.splitext(file)
    fullname = fn + ".agc"
    # fullname = fn + ".rms"
    print fullname

    # for i in range(0, len(d)):
    #     if abs(d[i]) > rms:
    #         m = i
    #         # d[i] = np.sign(d[i])*rms
    #
    # d[:m] = d[:m] / d[m] * rms
    # print m, d[m]
    #
    # tr.data = d

    tr.data = gain2(data=tr.data, npts=tr.stats.npts, delta=tr.stats.delta)

    tr.write(filename=fullname, format="SAC")

    # tr.data =


    # return rms

def time2depth_sac(file, dz, zrange, model="ak135"):

    print file

    tr = read(file)[0]
    dt = tr.stats.delta
    nt = tr.stats.npts

    t = np.arange(nt) * dt
    d = time2depth(time=t, type="P", model=model)

    # print t.shape, d.shape

    # print d[:10]
    # print tr.data.shape

    nz = int(round((zrange[1] - zrange[0]) / dz + 0.5))
    # print nz

    d_i = np.linspace(zrange[0], zrange[1], num=nz, endpoint=True)
    # print d_i.shape

    #  interpolate the depth profile into regular sampling interval
    f = interp1d(d, tr.data)
    s = f(d_i)

    tr.data = s
    tr.stats.delta = dz
    tr.stats.npts = nz

    tr.write(filename=file+".d", format="SAC")

    # os._exit(0)
    #
    # pass

def achalf(file):
    """
    cut the latter half part of the data


    :return:
    """

    tr = read(file)[0]
    npts = tr.stats.npts/2
    tr.data = tr.data[npts:]
    tr.stats.npts = len(tr.data)
    tr.write(file+"h", format="SAC")


    pass

##############################################################################
# plotting

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.transforms import blended_transform_factory
mpl.rcParams['axes.labelsize'] = 10
mpl.rcParams['axes.titlesize'] = 8
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['lines.markersize'] = 30
mpl.rcParams['lines.markeredgewidth'] = 0.8
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.size'] = 10

# plot the final section
# filelst: containning the station for each file in SAC format
# ofilename: output figure name
# tmin, tmax, the input file should be symetric,
# so that we cut from 180-250, that means 180 equals to zero .
# tzero: zero the first tzero time window to zero, aiming at reduce the large amplitude at t=0.
# if tzero < 0: then no zero action applied.
# agc: auto gain control, if agc > 0: then apply gain(), if agc < 0: then apply gain2(),
# if agc == 0: then no agc applied.
# taper: if taper < 0: no taper applied, if taper > 0: then taper is the window length.
# times = 1, apply times times taper
def plot_profile(filelst, ofilename, tmin=180, tmax=250, tzero=5, agc=5.0, taper=10.0, times=2.0, save=True, savepath="SDI/plot"):

    fig = plt.figure()

    filelst.sort()

    ll = len(filelst)

    # savepath = "SDI/plot"
    try:
        os.makedirs(savepath)
    except:
        pass

    # plot moho and LAB
    with open("lab.xyz", "r") as fp:
        lablst = fp.readlines()
    moho = []
    sta  = []
    lab = []
    for line in lablst:
        row = line.split()
        sta.append(row[1])
        moho.append(float(row[4]))
        lab.append(float(row[5]))

    tmoho = []
    tlab  = []
    for m in moho:
        tmoho.append(depth2time(m))
    for l in lab:
        tlab.append(depth2time(l))

    # print "length of tmoho:", len(tmoho), len(tlab)

    # read file
    st = read(filelst[0])
    tr = st[0]
    # tr.stats.sac.dist = 1
    tr.stats.distance = 1000
    # tr.filter("bandpass", freqmin=0.5, freqmax=4)

    t1 = tr.stats.starttime + tmin
    t2 = tr.stats.starttime + tmax
    tr.trim(starttime=t1, endtime=t2)
    tr.stats.starttime = 0

    # RMS and normalization
    # tr.normalize()
    # tr.data[np.where(tr.data>0.5)]  = 0.5
    # tr.data[np.where(tr.data<-0.5)] = -0.5
    # dd = tr.data
    # rms = math.sqrt(sum(dd*dd))/len(dd)
    # print rms
    # tr.data = tr.data/rms


    # zero the first tzero seconds
    if tzero>=0:
        n1 = 0
        n2 = int(round(float(tzero)/tr.stats.delta))
        tr.data[n1:n2] = 0.0

    # taper
    if taper > 0:
        hann = np.ones_like(tr.data)
        ntaper = int(round(taper / tr.stats.delta))
        hann1 = np.hanning(ntaper*2)
        hann[:ntaper] = hann1[:ntaper]
        hann = np.power(hann, times)
        tr.data *= hann
    else:
        pass

    # agc
    if agc<0:
        tr.data = gain2(tr.data, npts=tr.stats.npts, delta=tr.stats.delta)
    elif agc>0:
        tr.data = gain(tr.data, npts=tr.stats.npts, delta=tr.stats.delta, tlen=agc)
    else: # agc=0
        pass

    # save data
    bname = filelst[0]
    bname = savepath + "/" + os.path.basename(bname)

    stnm = tr.stats.sac.kstnm
    if stnm in sta:
        print stnm
        i = sta.index(stnm)
        tr.stats.sac.t5 = tr.stats.starttime + tmoho[i]
        tr.stats.sac.t6 = tr.stats.starttime + 0.5*(tmoho[i]+tlab[i])
        tr.stats.sac.t7 = tr.stats.starttime + tlab[i]

    tr.write(filename=bname, format="SAC")


    for i in range(1,len(filelst)):
        tr = read(filelst[i])[0]
        tr.stats.distance = (i+1)*1000

        # tr.filter("bandpass", freqmin=0.5, freqmax=4)

        t1 = tr.stats.starttime + tmin
        t2 = tr.stats.starttime + tmax
        tr.trim(starttime=t1, endtime=t2)
        tr.stats.starttime = 0.0

        if tzero>=0:
            n1 = 0
            n2 = int(round(float(tzero)/tr.stats.delta))
            tr.data[n1:n2] = 0.0

        # taper
        if taper > 0:
            hann = np.ones_like(tr.data)
            ntaper = int(round(taper / tr.stats.delta))
            hann1 = np.hanning(ntaper*2)
            hann[:ntaper] = hann1[:ntaper]
            hann = np.power(hann, times)
            tr.data *= hann
        else:
            pass

        if agc<0:
            tr.data = gain2(tr.data, npts=tr.stats.npts, delta=tr.stats.delta)
        elif agc>0:
            tr.data = gain(tr.data, npts=tr.stats.npts, delta=tr.stats.delta, tlen=agc)
        else: # agc = 0
            pass

        # save data
        bname = filelst[i]
        bname = savepath + "/" + os.path.basename(bname)

        stnm = tr.stats.sac.kstnm
        if stnm in sta:
            i = sta.index(stnm)
            tr.stats.sac.t5 = tr.stats.starttime + tmoho[i]
            tr.stats.sac.t6 = tr.stats.starttime + 0.5*(tmoho[i]+tlab[i])
            tr.stats.sac.t7 = tr.stats.starttime + tlab[i]

        tr.write(filename=bname, format="SAC")

        st += tr

    # print "reading file done."

    # st.filter("bandpass", freqmin=0.5, freqmax=4)

    st.plot(type="section", size=(1000,400), linewidth=0.5,
            time_down=True, fig=fig)

    plt.xlim(0, ll+2)
    plt.xlabel("")
    plt.xticks([])
    plt.yticks(range(0,75,5))

    # plot depth on the right side
    # ax = fig.axes[0]
    # a = [0, 14, 30, 49, 69, 89, 109, 130, 150, 170, 191, 211, 232, 253, 275]
    # cout = 0
    # transform = blended_transform_factory(ax.transData, ax.transAxes)
    # for i in range(0, 75, 5):
    #     ax.text(ll+1.5, i, str(a[cout]), va="center", ha="right")
    #     cout += 1
    # plt.text(ll+2.1, 35, "Depth [km]", rotation=270, va="center")

    a = [0, 14, 30, 49, 69, 89, 109, 130, 150, 170, 191, 211, 232, 253, 275]
    cout = 0
    for i in range(0, 75, 5):
        plt.text(ll+1.5, i, str(a[cout]), va="center", ha="right")
        cout += 1
    plt.text(ll+2.1, 35, "Depth [km]", rotation=270, va="center")

    # plot station name
    ax = fig.axes[0]
    transform = blended_transform_factory(ax.transData, ax.transAxes)
    for tr in st:
        ax.text(tr.stats.distance / 1e3, 1.0, tr.stats.station, rotation=270,
                va="bottom", ha="center", transform=transform, zorder=10)



    # the following commented lines moved forward so that write header into SAC
    # plot moho and LAB
    # with open("lab.xyz", "r") as fp:
    #     lablst = fp.readlines()
    # moho = []
    # sta  = []
    # lab = []
    # for line in lablst:
    #     row = line.split()
    #     sta.append(row[1])
    #     moho.append(float(row[4]))
    #     lab.append(float(row[5]))
    #
    # tmoho = []
    # tlab  = []
    # for m in moho:
    #     tmoho.append(depth2time(m))
    # for l in lab:
    #     tlab.append(depth2time(l))

    for tr in st:
        x = tr.stats.distance / 1e3
        stnm = tr.stats.sac.kstnm
        if stnm not in sta:
            continue

        i = sta.index(stnm)
        ax.plot(x, tmoho[i], "r_", markersize=20, alpha=0.5)
        try:
            ax.plot(x, tr.stats.t2, "m_", markersize=20, alpha=0.5)
        except:
            ax.plot(x, 0.5*(tmoho[i]+tlab[i]), "m_", markersize=20, alpha=0.5)
        ax.plot(x, tlab[i],  "c_", markersize=20, alpha=0.5)



    # fig.tight_layout()
    fig.savefig(ofilename)

    pass

# same as the plot_profile, but save the data in plot directory
# def plot_savedata(filelst, ofilename, tmin=180, tmax=250, tzero=5, agc=5.0, taper=10.0, times=2.0):
#     pass
def plot(filelst, ofilename, reverse=False):

    fig = plt.figure()

    # filelst.sort(reverse=reverse)

    ll = len(filelst)

    # read file
    st = read(filelst[0])
    tr = st[0]
    tr.stats.distance = 1000
    for i in range(1,len(filelst)):
        tr = read(filelst[i])[0]
        tr.stats.distance = (i+1)*1000
        st += tr

    # print "reading file done."

    # st.filter("bandpass", freqmin=0.5, freqmax=4)

    st.plot(type="section", size=(1500,150), linewidth=0.5,
            time_down=True, fig=fig)

    plt.xlim(0, ll+2)
    plt.xlabel("")
    plt.xticks([])
    plt.yticks(range(0,75,5))

    a = [0, 14, 30, 49, 69, 89, 109, 130, 150, 170, 191, 211, 232, 253, 275]
    cout = 0
    for i in range(0, 75, 5):
        plt.text(ll+1.5, i, str(a[cout]), va="center", ha="right")
        cout += 1
    plt.text(ll+2.1, 35, "Depth [km]", rotation=270, va="center")

    # plot station name
    ax = fig.axes[0]
    transform = blended_transform_factory(ax.transData, ax.transAxes)
    for tr in st:
        ax.text(tr.stats.distance / 1e3, 1.0, tr.stats.station, rotation=270,
                va="bottom", ha="center", transform=transform, zorder=10)

    for tr in st:
        x = tr.stats.distance / 1e3
        ax.plot(x, tr.stats.sac.t5, "r_", markersize=20, alpha=0.5)
        try:
            ax.plot(x, tr.stats.sac.t2, "m_", markersize=20, alpha=0.5)
        except:
            ax.plot(x, tr.stats.sac.t6, "m_", markersize=20, alpha=0.5)
        ax.plot(x, tr.stats.sac.t7, "c_", markersize=20, alpha=0.5)



    # fig.tight_layout()
    fig.savefig(ofilename)

    pass

def gain(data, npts, delta=0.025, tlen=5.0):

    nsec =  int(round(npts*delta/tlen))
    nt   = int(round(tlen/delta))

    # print nsec, nt

    x = []
    y = []
    xi = np.arange(len(data))
    # print min(xi), max(xi), len(xi)

    for i in range(nsec):
        # print i
        i1 = i*nt
        i2 = (i+1)*nt
        if i2 > npts:
            i2 = npts

        dd = data[i1:i2]
        rms = np.sqrt(sum(dd*dd))/len(dd)
        x.append(0.5*(i1+i2))
        y.append(rms)

        # print i1, i2, rms

    x1 = 0
    y1 = y[0]
    x.insert(0, x1)
    y.insert(0, y1)

    x2 = len(data)
    y2 = y[-1]
    x.append(x2)
    y.append(y2)

    x = np.array(x)
    y = np.array(y)

    # interpolation
    f = interp1d(x, y, kind="linear")
    g = f(xi)

    # print len(g), len(data)
    # print min(g), max(g)

    return data/g

    pass

# simple gain at fixed points
def gain2(data, npts, delta=0.025):

    tlen = (npts-1)*delta
    if tlen<70:
        print "data lenght is not enough. EXIT."
        os._exit(0)


    # i1 = 0
    # i2 = int(round(5.0/delta))
    # data[:i2] = 0.0

    # the first window 0-15 s, the fixed point at 7.5s
    tw1 = 0.0
    tw2 = 20.0

    x = []
    y = []

    i1 = int(round(tw1/delta))
    i2 = int(round(tw2/delta))

    dd = data[i1:i2]
    rms = np.sqrt(sum(dd*dd))/len(dd)

    x.append(0)
    y.append(rms)

    x.append(0.5*(i1+i2))
    y.append(rms)

    # the second window
    tw1 = 20.0
    tw2 = 40.0

    i1 = int(round(tw1/delta))
    i2 = int(round(tw2/delta))

    dd = data[i1:i2]
    rms = np.sqrt(sum(dd*dd))/len(dd)

    x.append(0.5*(i1+i2))
    y.append(rms)


    # the third window
    tw1 = 40.0
    tw2 = 70.0

    i1 = int(round(tw1/delta))
    i2 = int(round(tw2/delta))

    dd = data[i1:i2]
    rms = np.sqrt(sum(dd*dd))/len(dd)

    x.append(0.5*(i1+i2))
    y.append(rms)

    # the fourth window
    x.append(npts-1)
    y.append(rms)



    x = np.array(x)
    y = np.array(y)

    # print x
    # print y

    # interpolation
    xi = np.arange(len(data))
    f = interp1d(x, y, kind="linear")
    g = f(xi)

    # print min(g), max(g), npts, tlen

    # print len(g), len(data)
    # print min(g), max(g)

    return data/g



    # pass

def depth2time(depth, model="ak135", twoway=True):

    dep, vp, vs, tvp, tvs = d2t(model=model, twoway=twoway)
    f = interp1d(dep, tvp)

    return f(depth)

    pass

# save header that picks in plot folder
def export_header(filelst, file_header="header"):

    fp = open(file_header, "w")

    for file in filelst:

        # print file

        tr = read(file)[0]

        try:
            t5 = tr.stats.sac.t5
            t6 = tr.stats.sac.t6
            t7 = tr.stats.sac.t7
        except:
            continue

        try:
            t2 = tr.stats.sac.t2
            fp.write("%s\t%f\t%f\t%f\n"% (file, t5, t2, t7))
        except:
            t6 = tr.stats.sac.t6
            fp.write("%s\t%f\t%f\t%f\n"% (file, t5, t6, t7))


    fp.close()

# write header from files to sac files
def import_header(file_header="header", name=""):

    with open(file_header, "r") as fp:
        lst = fp.readlines()

    # print len(lst)

    for line in lst:

        # print line

        row = line.split()

        file = row[0]

        if name!="":
            file = file.replace("plot", "plot"+name)

        moho = float(row[1])
        mld  = float(row[2])
        lab  = float(row[3])

        try:
            tr = read(file)[0]
            tr.stats.sac.t5 = moho
            tr.stats.sac.t2 = mld
            tr.stats.sac.t7 = lab

            tr.write(file, format="SAC")
        except:
            pass



    pass


# end plotting
##############################################################################

if __name__=="__main__":

    # pre-processing
    # pathname = "processedSeismograms/Event_2007_10_26_09_17_04/*.sac"
    if True:
        pathname = "processedSeismograms/Event*/*BHZ.sac"
        lst = glob.glob(pathname)
        preprocess(filelst=lst, freqmin=0.05, freqmax=5.0, evdp_unit="m", sample_rate=40)

    # snr
    # pathname = "SDI/data/*/*BHZ.sac"
    # pathname = "SDI/data/1A.NE00/*.sac"
    if False:
        pathname = "SDI/data/*/*BHZ.sac"
        lst = glob.glob(pathname)
        snr(filelst=lst)

    # path, filename = filen(sacfile="/data/sun/work/nccsdi/data/1a/SDI/data/1A.NE00/Event_2007_10_26_09_17_04.BHN.sac")
    # print path
    # print filename

    # auto-correlation
    if False:
        pathname = "SDI/good/*/*BHZ.sac"
        lst = glob.glob(pathname)
        autocorrelate(filelst=lst)
