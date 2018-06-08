import numpy as np
from numpy import pi, polymul
from scipy import array, signal 
from scipy.signal import butter, lfilter, filtfilt
from collections.abc import Mapping
import time
import os, errno
import warnings


def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        return interpolator(x)

    def ufunclike(xs):
        return np.array(list(map(pointwise, array(xs))))

    return ufunclike


def stft(x, sr, fftsize=2 ** 10, overlap=4, ref=1, r_meas=1, r_new=1):
    hop = fftsize // overlap
    win = signal.hann(fftsize, sym=False) / (overlap // 2)
    padsize = int(fftsize * (np.ceil(len(x) / fftsize) - len(x) / fftsize))
    x = np.pad(x, (0, padsize), 'constant', constant_values=0)
    fftfreq = np.fft.rfftfreq(fftsize, 1 / sr)
    xfft = np.array([np.fft.rfft(win * x[i*hop:i*hop+fftsize])
                    for i in range(len(x) // hop - (fftsize // hop - 1))])
    xfft = np.abs(np.squeeze(xfft))

    r_meas = np.linalg.norm(r_meas-0)
    xdb = 20 * np.log10((r_meas/r_new)*xfft*2/fftsize/np.sum(win)/ref)
    return xfft, fftfreq, xdb, x


def istft(x, ifftsize=2 ** 10, overlap=4):
    hop = ifftsize // overlap
    win = signal.hann(ifftsize, sym=False) / (overlap // 2)
    xifft = np.zeros(np.size(x, 0) * hop + (overlap - 1) * hop)
    for n, i in enumerate(range(len(xifft) // hop - (ifftsize // hop - 1))):
        xifft[i * hop:i * hop + ifftsize] += win*(np.fft.irfft(x[n], ifftsize))  # overlap-add
    return xifft


def amp2db(x, ref=1):
    # return 20 * np.log10(np.clip(abs(x),2e-5, np.max(abs(x))) / ref)
    return 20 * np.log10(abs(x) / ref)


def db2amp(x, ref=1):
    return ref*10 ** (x / 20)


def dbadd(*args):
    db_sum = 0
    for i in args:
        db = np.asanyarray(i)
        db_sum += 10.0 ** (db / 10.0)
    db_sum = 10.0 * np.log10(db_sum)
    return db_sum

def dbavg(*args):
    db_avg = 0
    for i in args:
        db = np.asanyarray(i)
        db_avg += 10.0 ** (db / 10.0)
        db_avg /= len(args)
    db_avg = 10.0 * np.log10(db_avg)
    return db_avg


def dbsingle(dbspectrum):
    db_single = np.empty(np.size(dbspectrum, 0))
    for k in np.arange(np.size(dbspectrum, 0)):
        for i in np.arange(np.size(dbspectrum, 1)):
            db_single[k] += 10.0 ** (dbspectrum[k, i] / 10.0)
        db_single[k] = 10.0 * np.log10(db_single[k])
    return db_single


def a_weighting(x, fs):
    """
    Translated from a MATLAB script (which also includes C-weighting, octave
    and one-third-octave digital filters).
    Author: Christophe Couvreur, Faculte Polytechnique de Mons (Belgium)
        couvreur@thor.fpms.ac.be
    Last modification: Aug. 20, 1997, 10:00am.
    BSD license
    http://www.mathworks.com/matlabcentral/fileexchange/69
    Translated from adsgn.m to Python 2009-07-14 endolith@gmail.com
    """
    """Design of an A-weighting filter.
    b, a = A_weighting(fs) designs a digital A-weighting filter for
    sampling frequency `fs`. Usage: y = scipy.signal.lfilter(b, a, x).
    Warning: `fs` should normally be higher than 20 kHz. For example,
    fs = 48000 yields a class 1-compliant filter.
    References:
       [1] IEC/CD 1672: Electroacoustics-Sound Level Meters, Nov. 1996.
    """
    # Definition of analog A-weighting filter according to IEC/CD 1672.
    f1 = 20.598997
    f2 = 107.65265
    f3 = 737.86223
    f4 = 12194.217
    a1000 = 1.9997

    nums = [(2 * pi * f4) ** 2 * (10 ** (a1000 / 20)), 0, 0, 0, 0]
    dens = polymul([1, 4 * pi * f4, (2 * pi * f4) ** 2],
                   [1, 4 * pi * f1, (2 * pi * f1) ** 2])
    dens = polymul(polymul(dens, [1, 2 * pi * f3]),
                   [1, 2 * pi * f2])

    # Use the bilinear transformation to get the digital filter.
    # (Octave, MATLAB, and PyLab disagree about Fs vs 1/Fs)
    b, a = signal.bilinear(nums, dens, fs)
    y = signal.lfilter(b, a, x)
    return y


def a_weighting_f(f):
    f = np.array(f)
    f2 = np.power(f, 2)
    f4 = np.power(f, 4)
    ra = (12194. ** 2 * f4) / ((f2 + 20.6 ** 2) * np.sqrt((f2 + 107.7 ** 2) * (f2 + 737.9 ** 2)) * (f2 + 12194. ** 2))
    a = 20 * np.log10(ra) + 2
    return a


def rms_flat(a):  # from matplotlib.mlab
    """
    Return the root mean square of all the elements of *a*, flattened out.
    """
    return np.sqrt(np.mean(np.absolute(a) ** 2))


def spl_t(x, sr, interval=0.125, pref=2e-05):
    pref = rms_flat(pref)
    n = int(sr * interval)
    padsize = int(n * (np.ceil(len(x) / n) - len(x) / n))
    x = np.pad(x, (0, padsize), 'constant', constant_values=0)
    level = np.empty(int(len(x) / n))
    level_a = np.empty(int(len(x) / n))
    t = np.arange(int(len(x) / n)) * interval
    x_a = a_weighting(x, sr)
    for i in range(int(len(x) / n)):
        level[i] = 20 * np.log10(rms_flat(x[i * n:(i + 1) * n]) / pref)
        level_a[i] = 20 * np.log10(rms_flat(x_a[i * n:(i + 1) * n]) / pref)
    return level, level_a, t


def lpf(x, fs, tau=0.125):
    alpha = 1 / (tau * fs)  # approximation of alpha = 1 - np.exp(-1 / (tau*fs))
    y = np.zeros_like(x)
    xa = a_weighting(x, fs)
    ya = np.zeros_like(xa)
    yk = x[0]
    yka = xa[0]
    for k in range(len(y)):
        yk += alpha * (x[k] - yk)
        y[k] = yk
        yka += alpha * (xa[k] - yka)
        ya[k] = yka
    return y, ya



def fspace(fc=np.array([63., 125., 250., 500., 1000., 2000., 4000., 8000.]),
           nf=15, bw=1, sc="log"):
    nb = len(fc)  # number of bands
    flimits = np.zeros((nb, 2))  # band frequency limits
    if bw == 1:
        for k in range(nb):
            flimits[k, 0] = (fc[k] / (np.sqrt(2))) + 1
            flimits[k, 1] = (fc[k] * (np.sqrt(2))) - 1
    elif bw == 3:
        for k in range(nb):
            flimits[k, 0] = fc[k] / (2 ** (1 / 6))
            flimits[k, 1] = fc[k] * (2 ** (1 / 6))
    else:
        Exception('Non-supported argument for input "bandWidth"')
    freqvec = np.zeros((nb, nf))
    if sc == 'log':
        for k in range(nb):
            freqvec[k] = np.logspace(np.log10(flimits[k, 0]),
                                     np.log10(flimits[k, 1]), nf)
    elif sc == 'lin':
        for k in range(nb):
            freqvec[k] = np.linspace(flimits[k, 0], flimits[k, 1], nf)
    else:
        Exception('Non-supported argument for input "scale"')
    return freqvec, flimits


def dbsum(levels, axis=None):
    """Energetic summation of levels.

    :param levels: Sequence of levels.
    :param axis: Axis over which to perform the operation.

    .. math:: L_{sum} = 10 \\log_{10}{\\sum_{i=0}^n{10^{L/10}}}

    """
    levels = np.asanyarray(levels)
    return 10.0 * np.log10((10.0**(levels/10.0)).sum(axis=axis))


def narrow2third(xdb, fnarrow):
    fc_pref = np.array([12.5, 16., 20., 25., 31.5, 40., 50., 63., 80., 100., 125.,
                        160., 200., 250., 315., 400., 500., 630., 800., 1000.,
                        1250., 1600., 2000., 2500., 3150., 4000., 5000., 6300.,
                        8000., 10000., 12500., 16000., 20000.])
    # determine lower and upper limits of each 1/3 octave band
    fc = np.zeros(len(fc_pref))
    fc_lim = np.zeros((2, len(fc_pref)))
    bands = np.zeros(len(fc_pref))
    for a in range(len(fc_pref)):
        fc[a] = 1000*(2**(1/3))**(a-19)
        fc_lim[0, a] = fc[a]/2**(1/6)
        fc_lim[1, a] = fc[a]*2**(1/6)
        idx = np.where(np.logical_and(
            fnarrow >= fc_lim[0, a], fnarrow < fc_lim[1, a]))
        idx = idx[0]
        bands[a] = 0
        if np.size(idx) == 0:
            warnings.warn('no point found in band centered at %i' % fc_pref[a])
        elif np.size(idx) == 1:
            warnings.warn('only one point found in band centered at %i' % fc_pref[a])
            bands[a] = xdb[idx]
        else:
            bands[a] = dbsum(xdb[idx])
    return bands, fc_pref
