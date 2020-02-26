import numpy as np
import math
from scipy.optimize import curve_fit
from scipy.signal import sosfilt, butter
import scipy.signal as spSig
import scipy as sp
import scipy.stats


def amp2db(x, ref=1):
    return 20 * np.log10(abs(x) / ref)


def NormAmp2dB(data):
    '''Converts amplitude to dB, using the maximum
    value as reference'''
    return 10 * np.log10(abs((data)/max(data)))


def find_nearest(array, value):
    '''
    Compares the input array with a the input value and outputs the index
    that corresponds to the nearest value of the input array
    '''
    idx = (np.abs(array-value)).argsort()[:2]
    return idx


def dB_clipp(data, threshold):
    '''
    Cuts a signal below a certain threshold
    '''
    Level = 20 * np.log10(abs(data))
    cut = np.where(Level > threshold)
    out = [Level[cut], data[cut], cut]

    return out


def frequencyBand(x, b=1, G=10, fref=1000.0):
    """
    Calculation of exact center frequeny and bandedge frequencies
    acording to IEC 61260:1995

    x: band index. typical ranges:
                                octave bands = 24 to 34
                                third octave bands 11 to 43
    b: see b in standard (inverse of bandwidth designator)
    G: octave ratio, only base 10 and 2 are supported
    fref: reference frequency, usualy 1000 Hz

    returns fm - midband frequency
            f1 - lower bandedge frequency
            f2 - upper bandedge frequency
    """
    if G == 10:
        G = 10 ** (3 / 10)
    if b % 2 == 1:  # odd
        fm = fref * G**((x-30) / (b))
    else:               # even
        fm = fref * G**((2 * x - 59) / (2 * b))

    f1 = fm * (G ** (-1 / (2 * b)))
    f2 = fm * (G ** (1 / (2 * b)))
    return fm, f1, f2


def round_sigfigs(num, sig_figs):
    """
    Round to specified number of sigfigs.
    """
    if num != 0:
        return round(
            num,
            -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0  # Can't take the log of 0


def lorentz(x, *p):
    '''
    Creates a Lorentzian curve at the specified x vector given the values:
    I - magnitude of the center frequency
    gamma - the gamma value of the curve
    x0 - the center frequency
    '''
    I, gamma, x0 = p

    return I * gamma**2 / ((x - x0)**2 + gamma**2)


def fit_lorentz(p, x, y):
    '''
    Fits a Lorentzian curve created using the Lorentz function above,
    to  the given data x,y
    '''
    return curve_fit(lorentz, x, y, p0=p)


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h


def half_hann(length, side='left'):
    '''Creates a half Hanning window
    Arguments: length - the length of the window in samples
    side - the side of the window's fade out. values: left or right
    '''
    win_samples = int(length)
    tvec = np.linspace(0, np.pi/2, win_samples)
    if side == 'right':
        win = np.cos(tvec)**2
    elif side == 'left':
        win = np.sin(tvec)**2
    return win


def butterWorthFilter(in_data,
                      sr,
                      frequencies=[],
                      order=3,
                      filter_type='highpass',
                      method='ba'):
    '''Creates a Butterworth filter and filters the given data
    Inputs:
    - in_data: input signal.
    - sr: sample rate.
    - frequencies: cutoff frequency for the highpass or lowpass filters or
    [fmin, fmax] for the bandpass and bandstop.
    - order: filter order.
    - filter_type: 'highpass', 'lowpass', 'bandpass' or bandstop'.
    - method: 'ba', zpk' or 'sos'. This corresponds to the creation of the
    filter. The filter is always converted to second order segments prior to
    filtering.

    Output:
    -out_data: the filtered output signal
    '''

    supported_filters = ['highpass', 'lowpass', 'bandpass', 'bandstop']
    supported_methods = ["ba", "zpk", "sos"]
    if not isinstance(frequencies, list):
        frequencies = [frequencies]

    # Input error check
    if filter_type not in supported_filters:
        raise ValueError('Unsupported filter type. Available'
                         + str(supported_filters))
    if any([isinstance(ck, str) for ck in frequencies]):
        raise ValueError("Frequency/ies must be int of float")
    if filter_type in supported_filters[:2] and len(frequencies) != 1:
        raise ValueError(filter_type + " filter requires only one fc.")
    elif filter_type in supported_filters[2:] and len(frequencies) != 2:
        raise ValueError(filter_type
                         + " filter requires a list of critical"
                         + " frequencies given as [f0, f1].")
    if method not in supported_methods:
        raise ValueError("Unknown method for filter design. "
                         + "Please select from: " + str(supported_methods))

    # Filter design
    nyq = 0.5 * sr
    if len(frequencies) == 1:
        fc = frequencies[0] / nyq
        w = np.array(fc)
    elif len(frequencies) == 2:
        low = frequencies[0] / nyq
        high = frequencies[1] / nyq
        w = np.array([low, high])
    if method == 'ba':
        b, a = butter(order, w, btype=filter_type, output='ba')
        sos = spSig.tf2sos(b, a)
    elif method == 'zpk':
        z, p, k = butter(order, w, btype=filter_type, output='zpk')
        sos = spSig.zpk2sos(z, p, k)
    elif method == 'sos':
        sos = butter(order, w, btype=filter_type, output='sos')
    out_data = sosfilt(sos, in_data)

    return out_data


def band_filter(data, sr, f_min, f_max, filt_order=3, bandWidth="third"):
    '''Octave/Third-octave band filter.

    Arguments: f_min,f_max - min/max frequency band of interest. Values are
                              truncated to the closest center frequency.
                bandWidth - chooses between octave and 3rd octave bands

    Returns a list with the  band filtered data and the corresponding
    frequency vector
    '''

    centerFreqs_3rd = np.array([12.5, 16.0, 20.0, 25.0, 31.5, 40.0, 50.0, 63.0,
                                80.0, 100.0, 125.0, 160.0, 200.0, 250.0, 315.0,
                                400.0, 500.0, 630.0, 800.0, 1000.0, 1250.0,
                                1600.0, 2000.0, 2500.0, 3150.0, 4000.0,
                                5000.0, 6300.0, 8000.0, 10000.0, 12500.0,
                                16000.0, 20000.0])

    centerFreqs_oct = np.array([16.0, 31.5, 63.0, 125.0, 250.0, 500.0,
                                1000.0, 2000.0, 4000.0, 8000.0, 16000.0])

    output = []
    if f_min > f_max:
        raise RuntimeError("Minimum frequency greater than maximum frequency")
    if bandWidth == "third":
        centerFreqs = centerFreqs_3rd
        ISO_f = np.array([frequencyBand(fc, 3) for fc in range(11, 44)])
    elif bandWidth == "one":
        centerFreqs = centerFreqs_oct
        ISO_f = np.array([frequencyBand(fc, 1) for fc in range(24, 34)])
    else:
        raise RuntimeError("Unsupported frequency bandWidth. "
                           + " Current options: \"third\", \"one\" ")
    f_min = min(ISO_f[:, 0], key=lambda x: abs(x-f_min))
    f_max = min(ISO_f[:, 0], key=lambda x: abs(x-f_max))
    f_min_idx = min(np.where(ISO_f[:, 0] == f_min)).item()
    f_max_idx = min(np.where(ISO_f[:, 0] == f_max)).item()
    centerFreqs_out = centerFreqs[f_min_idx: f_max_idx + 1]
    ISO_f_out = ISO_f[f_min_idx:f_max_idx+1, :]
    for n in range(len(ISO_f_out[:, 0])):
        f_low = ISO_f_out[n, 1]
        f_high = ISO_f_out[n, 2]
        filtered_data = butterWorthFilter(
            data,
            sr,
            frequencies=[f_low, f_high],
            order=filt_order,
            filter_type='bandpass',
            method='sos'
        )
        output.append(filtered_data)
    return output, centerFreqs_out
