import numpy as np
from numpy import pi, polymul
from scipy import array, signal
from collections.abc import Mapping
import time
# from scipy.interpolate import interp1d
import warnings

def sftpget(remotefile, remotepath='/home/gz/peaks_calc/', localpath='',
            adr="tapc18.ta.chalmers.se", usr="gz", pwd="zaxosZAXOS0"):
    import paramiko
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(adr, username=usr, password=pwd)
    sftp = ssh.open_sftp()
    localpath = localpath + remotefile
    remotepath = remotepath + remotefile
    sftp.get(remotepath, localpath)
    sftp.close()
    ssh.close()


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


def cnossos_source(speed):
    fc = np.array([63., 125., 250., 500., 1000., 2000., 4000., 8000.])
    ref_speed = 70
    coeffs = np.zeros([3, 6, 8])
    # category 1 (passenger cars)
    coeffs[0, :, :] = [[79.7, 85.7, 84.5, 90.2, 97.3, 93.9, 84.1, 74.3],  # A_R
                       [30.0, 41.5, 38.9, 25.7, 32.5, 37.2, 39.0, 40.0],  # B_R
                       [94.5, 89.2, 88.0, 85.9, 84.2, 86.9, 83.3, 76.1],  # A_P
                       [-1.3, 7.2, 7.7, 8.0, 8.0, 8.0, 8.0, 8.0],  # B_P  # b
                       [0.0, 0.0, 0.0, 2.6, 2.9, 1.5, 2.3, 9.2],  # a
                       [0.0, 0.0, 0.0, -3.1, -6.4, -14.0, -22.4, -11.4]]  # b

    # category 2 (medium heavy vehicles)
    coeffs[1, :, :] = [[84.0, 88.7, 91.5, 96.7, 97.4, 90.0, 83.8, 80.5],  # A_R
                       [30.0, 35.8, 32.6, 23.8, 30.1, 36.2, 38.3, 40.1],  # B_R
                       [101., 96.5, 98.8, 96.8, 98.6, 95.2, 88.8, 82.7],  # A_P
                       [-1.9, 4.7, 6.4, 6.5, 6.5, 6.5, 6.5, 6.5],  # B_P
                       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # a
                       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]  # b

    # category 3 (heavy duty vehicles)
    coeffs[2, :, :] = [[87.0, 91.7, 94.1, 100.7, 100.8, 94.3, 87.1, 82.5],  # A_R
                       [30.0, 33.5, 31.3, 25.4, 31.8, 37.1, 38.6, 40.6],  # B_R
                       [104.4, 100.6, 101.7, 101.0, 100.1, 95.9, 91.3, 85.3],  # A_P
                       [0.0, 3.0, 4.6, 5.0, 5.0, 5.0, 5.0, 5.0],  # B_P
                       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # a
                       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]  # b

    lw_roll = coeffs[:, 0, :] + coeffs[:, 1, :] * np.log10(speed / ref_speed)
    lw_prop = coeffs[:, 2, :] + coeffs[:, 3, :] * ((speed - ref_speed) / ref_speed)
    lw_ex = np.zeros([3, 11])

    # for i in range(3):
    #     lw_roll_ex[i, :] = np.insert(np.array(lw_roll[i, :]), [0, 0, len(lw_roll[i, :])], interp([16., 31.5, 16000.], fc, lw_roll[i, :]))
    #     lw_prop_ex[i, :] = np.insert(np.array(lw_prop[i, :]), [0, 0, len(lw_roll[i, :])], interp([16., 31.5, 16000.], fc, lw_prop[i, :]))
    fc = np.array([16., 31.5, 63., 125., 250., 500., 1000., 2000., 4000., 8000., 16000.])

    lw = dbadd(lw_roll, lw_prop)

    t = 90
    a = (lw[..., 0]-t)/fc[2]
    b = fc[:2]
    lw_ex[..., :2] = np.dot(a[..., None], b[None, ...])+t
    lw_ex[..., 2:-1] = lw
    lw_ex[..., -1] = lw[..., -1]
    return lw_ex, fc


def imagine_source_new(speed):
    ref_speed = 70
    coeffs = np.zeros([3, 4, 27])
    # category 1 (passenger cars)
    coeffs[0, :, :] = [[69.9, 69.9, 69.9, 74.9, 74.9, 74.9, 79.3, 82, 81.2, 80.9, 78.9, 78.8, 80.5, 85, 87.9, 90.9, 93.3, 92.8, 91.5, 88.5, 84.9, 81.8, 78.7, 74.9, 71.8, 69.1, 65.6],  # A_R
                       [33, 33, 33, 30, 30, 30, 41, 41.2, 42.3, 41.8, 38.6, 35.5, 32.9, 25, 25, 27, 33.4, 36.7, 37, 37.5, 37.5, 38.6, 39.6, 40, 39.9, 40.2, 40.3],  # B_R
                       [87, 87, 87, 87.9, 90.8, 89.9, 86.9, 82.6, 81.9, 82.3, 83.9, 83.3, 82.4, 80.6, 80.2, 77.8, 78, 81.4, 82.3, 82.6, 81.5, 80.2, 78.5, 75.6, 73.3, 71, 68.1],  # A_P
                       [0, 0, 0, 0, -3, 0, 8, 6, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]]  # B_P

    # category 2 (medium heavy vehicles)
    coeffs[1, :, :] = [[76.5, 76.5, 76.5, 78.5, 79.5, 79.5, 82.5, 84.3, 84.7, 84.3, 87.4, 87.8, 89.8, 91.6, 93.5, 94.6, 92.4, 89.6, 88.1, 85.9, 82.7, 80.7, 78.8, 76.8, 76.7, 75.7, 74.5],
                       [33, 33, 33, 30, 30, 30, 32.9, 35.9, 38.1, 36.5, 33.5, 30.6, 27.7, 21.9, 23.8, 28.4, 31.1, 35.4, 35.9, 36.7, 36.3, 37.7, 38.5, 39.8, 39.9, 40.2, 40.3],
                       [93.9, 93.9, 94.1, 95, 97.3, 96.1, 92.5, 91.9, 90.4, 93.4, 94.4, 94.2, 93, 90.8, 92.1, 92.5, 94.1, 94.5, 92.4, 90.1, 87.6, 85.8, 83.8, 81.4, 80, 77.2, 75.4],
                       [0, 0, 0, 0, -4, 0, 4, 5, 5.5, 6, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5]]

    # category 3 (heavy duty vehicles)
    coeffs[2, :, :] = [[79.5, 79.5, 79.5, 81.5, 82.5, 82.5, 85.5, 87.3, 87.7, 87.3, 89.5, 90.5, 93.8, 95.9, 97.3, 98, 95.6, 93.2, 91.9, 88.9, 85.5, 84.1, 82.2, 79.8, 78.6, 77.5, 76.8],
                       [33, 33, 33, 30, 30, 30, 31.4, 32.8, 36, 34.6, 32.7, 29.3, 26.4, 24.2, 25.9, 30.4, 32.3, 36.5, 36.8, 38, 36.8, 38.5, 38.9, 38.5, 40.2, 40.8, 41],
                       [95.7, 94.9, 94.1, 96.8, 101.8, 98.6, 95.5, 96.2, 95.7, 97.2, 96.3, 97.2, 95.8, 95.9, 96.8, 95.1, 95.8, 95, 92.7, 91.2, 88.7, 87.6, 87.2, 84.2, 82.7, 79.7, 77.6],
                       [0, 0, 0, -4, 0, 4, 3, 3, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]]

    n_axles = np.array([[5], [5], [6]])
    # dL_DAC16=0.95,  % addition due to asphalt DAC (Dense Asphalt Concrete), 0/16 (Stone sizes from 0 to 16 mm), 2 years old
    dl_dac16 = 0
    # lw_roll = dl_dac16 + 10*np.log10(n_axles/5)
    lw_roll =  coeffs[:, 0, :] + coeffs[:, 1, :] * np.log10(speed / ref_speed) + 10*np.log10(n_axles/5)
    lw_prop = coeffs[:, 2, :] + coeffs[:, 3, :]*((speed - ref_speed)/ref_speed)
    fc = np.array([25., 31.5, 40., 50., 63., 80., 100., 125., 160.,
                   200., 250., 315., 400., 500., 630., 800., 1000.,
                   1250., 1600., 2000., 2500., 3150., 4000., 5000.,
                   6300., 8000., 10000.])
    lw_ex = np.zeros([3, 33])
    lw_roll_ex = np.zeros([3, 33])
    lw_prop_ex = np.zeros([3, 33])

    # for i in range(3):
    #     lw_roll_ex[i, :] = np.insert(np.array(lw_roll[i, :]), [0, 0, 0, len(lw_roll[i, :]), len(lw_roll[i, :]), len(lw_roll[i, :])], interp([12.5, 16., 20., 12500., 16000., 20000.], fc, lw_roll[i, :]))
    #     lw_prop_ex[i, :] = np.insert(np.array(lw_prop[i, :]), [0, 0, 0, len(lw_roll[i, :]), len(lw_roll[i, :]), len(lw_roll[i, :])], interp([12.5, 16., 20., 12500., 16000., 20000.], fc, lw_prop[i, :]))
    #     # fi_roll = interp1d(fc, lw_roll[i, :])
    #     # fi_prop = interp1d(fc, lw_prop[i, :])
    #     # fx_roll = extrap1d(fi_roll)
    #     # fx_prop = extrap1d(fi_prop)
    #     # lw_roll_ex[i,:] = np.insert(np.array(lw_roll[i, :]), 0, fx_roll([12.5, 16., 20.]))
    #     # lw_prop_ex[i,:] = np.insert(np.array(lw_prop[i, :]), 0, fx_prop([12.5, 16., 20.]))

    lw = dbadd(lw_roll, lw_prop)
    fc = np.array([12.5, 16., 20., 25., 31.5, 40., 50., 63., 80., 100., 125.,
                   160., 200., 250., 315., 400., 500., 630., 800., 1000.,
                   1250., 1600., 2000., 2500., 3150., 4000., 5000., 6300.,
                   8000., 10000., 12500., 16000., 20000.])

    b = fc[:3]

    t = 40
    a = (lw_roll[..., 0]-t)/fc[3]
    lw_roll_ex[..., :3] = np.dot(a[..., None], b[None, ...])+t
    lw_roll_ex[..., 3:-3] = lw_roll
    lw_roll_ex[..., -3:] = np.tile(lw_roll[..., -1:], 3)

    t = 80
    a = (lw_prop[..., 0]-t)/fc[3]
    lw_prop_ex[..., :3] = np.dot(a[..., None], b[None, ...])+t
    lw_prop_ex[..., 3:-3] = lw_prop
    lw_prop_ex[..., -3:] = np.tile(lw_prop[..., -1:], 3)

    lw_ex = dbadd(lw_roll_ex, lw_prop_ex)
    return lw_ex, lw_prop_ex, lw_roll_ex, fc


def airattdb1m(fvec, ht=40, tempC=24, pa=101325):
    pr = 101325  # Reference ambient pressure
    t = 273.15 + tempC  # Temperature in Kelvin
    t0 = 293.15  # Reference temperature
    t01 = 273.16

    c = -6.8346 * (t01 / t) ** 1.261 + 4.6151
    h = ht * pr / pa * 10 ** c  # Molar conc.of water vapour( %)
    # Relaxation freq of oxygen
    fro = (pa / pr) * (24 + 4.04 * 10 ** 4 * h * (0.02 + h) / (0.391 + h))
    frn = (pa / pr) / np.sqrt(t / t0) * (
        9 + 280 * h * np.exp(-4.17 * ((t / t0) ** (-1 / 3) - 1)))
    airatt = 8.686 * fvec ** 2. * (1.84 * 1e-11 * (pr / pa) * np.sqrt(t / t0)
                                   + (t / t0) ** (-2.5)
                                   * (0.01275 * np.exp(-2239.1 / t)
                                   * (fro + fvec ** 2 / fro) ** (-1)
                                   + 0.1068 * np.exp(-3352 / t)
                                   * (frn + fvec ** 2 / frn) ** (-1)))
    return airatt


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


def stddb(veh, speed, ver='new'):
    speed = np.array(speed)
    if ver == 'new':
        if veh == 0:
            fac1 = 4.94
            fac2 = -0.7
        if veh == 1:
            fac1 = 6.35
            fac2 = -0.9
            speed[speed < 50] = 50
    if ver == 'old':
        if veh == 0:
            fac1 = 5.5
            fac2 = -0.7
        if veh == 1:
            fac1 = 10
            fac2 = -0.9
            speed[speed < 50] = 50
    y = fac1*np.exp(fac2*speed/50)
    return (np.random.random()-(y/2))*2


def track_job(job, update_interval=1):
    while job._number_left > 0:
        print("{0}".format(job._number_left * job._chunksize),
              end='|', flush=True)
        time.sleep(update_interval)


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


def rec_key_replace(obj):
    if isinstance(obj, Mapping):
        return {key.replace('/', '_'): rec_key_replace(val) for key, val in obj.items()}
    return obj

def saveToMat(obj):
    # obj = obj.item()
    mat = rec_key_replace(obj)
    return mat
