import numpy as np
import utils as gz
from scipy.optimize import curve_fit
from scipy.signal import butter, lfilter, filtfilt, sosfiltfilt, sosfilt, freqz, sosfreqz
import scipy.signal as spSig
import scipy as sp
import scipy.stats


def NormAmp2dB(data):
    '''Converts amplitude to dB, using the maximum
    value as reference'''
    return 10 * np.log10(abs((data)/max(data)))

def removeDelay(data, sr):
    '''
    Removes the initial delay introduced in an Impulse response, due to the source receiver distance
    or other delays introduced in the measurement chain.
    '''
    IR_a = gz.a_weighting(data, sr)
    xcross_idx = np.where(np.diff(np.sign(IR_a)))[0]
    avg = np.empty(len(xcross_idx)-1)

    for i, _ in enumerate(xcross_idx[1:]):
        avg[i] = np.average(abs(data[xcross_idx[i]:xcross_idx[i+1]+1]))

    start_idx = xcross_idx[np.argmax(avg > 2e-4)]

    IR_head = data[:start_idx]
    IR_tail = data[start_idx:]

    return IR_head, IR_tail, start_idx


def find_nearest(array,value):
    '''
    Compares the input array with a the input value and outputs the index that corresponds to the nearest
    value of of the input array
    '''
    idx = (np.abs(array-value)).argsort()[:2]
    return idx


def dB_clipp(data, threshold):
    '''
    Cuts a signal below a certain threshold
    '''
    Level = 20 * np.log10(abs(data))
    cut = np.where(Level > thershold)
    out = [Level[cut],data[cut],cut]

    return out

def findTau(data, tVec, sr, fit_start, fit_stop, threshold):
    '''Finds the end of exponential decay (referenced as the time "tau"
    in the Schroeder curve), setting the correct limit for the backwards
    integration method.

    The method creates a straight line that scans the data signal and
    calculates the rms difference between the signal and each line.
    The lowest rms value will correspond to the point of the curve where
    the exponential decay has ended.

    It is assumed that the data are sufficiently smooth (using band filtering,
    Hilbert transform and moving average filters) and converted to dB using the
    method amp2dB.

    It is also assumed that the data consist of two parts. An exponential decay
    and a steady background noise level.

    The fitted straight line is created form the following points:
    maxPoint: the maximum value of the data if the y axis
    threshold: the value of the data in the y axis that corresponds to the
               expected noise floor
    tVec[maxPoint]: the x value of the maxPoint
    tVec[testPoint]: the point that moves along the x axis to scan the data

    The filst 2000 values of the data are not scanned since the beggining of
    the impulse responce consists of zeros, due to the source - receiver distance
    in the measurement procedure. A SCRIPT removing this will be implemented in the
    future.
    '''

    maxPoint = np.argmax(data)
    conVal = []
    rmsdiff = []
    testPoints = range(fit_start, fit_stop, int(0.001 * sr))
    for testPoint in testPoints:
        b = (threshold - data[maxPoint])/(tVec[testPoint]-tVec[maxPoint])
        y2 = b*tVec[maxPoint:testPoint] + data[maxPoint]
        y = data[maxPoint:testPoint]

        rmsdiff.append(np.sqrt(((y - y2) ** 2).sum() / len(y)))
    tau = testPoints[np.argmin(rmsdiff)]

    return tau

def ampDistribution(data, tVec):
    values = np.array(np.sort(data))
    tau = findTau(values, tVec, values[0])
    noise_amplitude = values[tau]

    return noise_amplitude

def smooth(x,window_len=11,window='hanning'):
        """smooth the data using a window with requested size.

        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.

        input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

        output:
        the smoothed signal

        NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

        if x.ndim != 1:
            raise ValueError("smooth only accepts 1 dimension arrays.")
        if x.size < window_len:
            raise ValueError("Input vector needs to be bigger than window size.")
        if window_len<3:
            return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

        s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
        if window == 'flat': #moving average
            w=np.ones(window_len,'d')
        else:
            w=eval('np.'+window+'(window_len)')

        y=np.convolve(w/w.sum(),s,mode='valid')

        return y[int(window_len/2-1):-int(window_len/2)]
        # return y

def lorentz(x, *p):
    '''Creates a Lorentz destribution using I, gamma and x0 values'''
    I, gamma, x0 = p

    return I * gamma**2 / ((x - x0)**2 + gamma**2)

def fit_lorentz(p, x, y):
    '''Finds the coeficients for the Lorentz distribution that
        fits best to the data x,y'''
    return curve_fit(lorentz, x, y, p0 = p)

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

def half_hann(length, side='left'):
    '''Creates a half hanning window
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



def butterWorthFilter(in_data, sr, frequecies=[], order=3, filter_type='highpass', method='ba'):
    '''Creates a Butterworth filter and filters the given data
    Inputs:
    - in_data: input signal.
    - sr: sample rate.
    - frequencies: cutoff frequency for the highpass or lowpass filters or [fmin, fmax] for the bandpass and bandstop.
    - order: filter order.
    - filter_type: 'highpass', 'lowpass', 'bandpass' or bandstop'.
    - method: 'ba', zpk' or 'sos'. This corresponds to the creation of the filter. The filter is always
    converted to second order segments prior to filtering.

    Output:
    -out_data: the filtered output signal
    '''

    if filter_type != 'highpass' and filter_type != 'lowpass' and filter_type != 'bandpass' and filter_type != 'bandstop':
        print("wrong filter type. Please define filter type from the following options:")
        print("highpass, lowpass, bandstop, bandpass")
    nyq = 0.5 * sr
    if len(frequecies) == 1:
        fc = frequecies[0] / nyq
        w = np.array(fc)
    elif len(frequecies) == 2:
        low = frequecies[0] / nyq
        high = frequecies[1] / nyq
        w = np.array([low, high])
    else:
        print("wrong frequency input. Please define the frequecies list as:")
        print("[f_low, f_high], for bandpass and bandstop filters,")
        print("[f_c] for highpass and lowpass filters")
    if method == 'ba':
        b, a = butter(order, w, btype=filter_type,output='ba' )
        sos = spSig.tf2sos(b, a)
    elif method == 'zpk':
        z, p, k = butter(order, w, btype=filter_type, output='zpk')
        sos = spSig.zpk2sos(z, p,k)
    elif method == 'sos':
        sos = butter(order, w, btype=filter_type, output='sos')
    else:
        print("Unknown method for filter desigh. Please define the method parameter as:")
        print("ba or zpk or sos")
    out_data = sosfilt(sos, in_data)

    return out_data

def band_filter(data, sr, f_min, f_max, filt_order=3, bandWidth = "third"):
    '''Octave/Third-octave band filter.

    Arguements: f_min,f_max - min/max frequency band of interest. Values are
                              turncated to the closest center frequency.
                bandWidth - chooses between octave and 3rd octave bands

    Returns a list with the  band filtered data and the corresponding frequency vector
    '''

    centerFreqs_3rd = np.array([12.5, 16.0, 20.0, 25.0, 31.5, 40.0, 50.0, 63.0,
                                80.0, 100.0, 125.0, 160.0, 200.0, 250.0, 315.0,
                                400.0,500.0, 630.0, 800.0, 1000.0, 1250.0, 1600.0,
                                2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0,
                                8000.0, 10000.0, 12500.0, 16000.0, 20000.0])

    centerFreqs_oct = np.array([16.0, 31.5, 63.0, 125.0, 250.0, 500.0,
                                1000.0, 2000.0, 4000.0, 8000.0, 16000.0])


    output = []
    if f_min > f_max: raise RuntimeError("Minimum frequency greater than maximum frequency")
    if bandWidth == "third":
        centerFreqs = centerFreqs_3rd
        band_factor = 3
    elif bandWidth == "one":
        centerFreqs = centerFreqs_oct
        band_factor = 1
    else:
        raise RuntimeError("Unsupported frequency bandWidth. Current options: \"third\", \"one\" ")
    f_min = min(centerFreqs, key=lambda x:abs(x-f_min))
    f_max = min(centerFreqs, key=lambda x:abs(x-f_max))
    f_min_idx = min(np.where(centerFreqs == f_min)).item()
    f_max_idx = min(np.where(centerFreqs == f_max)).item()
    centerFreqs_out = centerFreqs[f_min_idx:f_max_idx+1]
    for f_c in centerFreqs_out:
        f_low = f_c * 2 ** (-1 / ( 2 * band_factor ))
        f_high = f_c * 2 ** (1 / ( 2 * band_factor ))
        filtered_data = butterWorthFilter(data, sr, frequecies=[f_low, f_high], order=filt_order, filter_type='bandpass', method='sos')
        output.append(filtered_data)
        centerFreqs = centerFreqs

    return output, centerFreqs_out


