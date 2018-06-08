import numpy as np
import utils as gz
from scipy.optimize import curve_fit
from scipy.signal import butter, lfilter, filtfilt, sosfiltfilt, sosfilt

class filters(object):

    ''' Various processing for time signals
    Properties - data: time signal
               - sr: the sample rate
               - output: list used for output


    Methods
    band_filter(f_min, f_max, bandWidth): an Octave/Third-octave band filter.

    Arguements: f_min,f_max - min/max frequency band of interest. Values are
                              turncated to the closest center frequency.
                bandWidth - chooses between octave and 3rd octave bands

    Returns a list with the filtered data

    butter_bandpass(f_low, f_high): returns the filter coeficients for a second
                                    order bandpass butterworth filtered


    butter_bandpass_filter(f_low, f_high): filters the signal using the filter
                                           created above
    '''
    centerFreqs_3rd = np.array([12.5, 16.0, 20.0, 25.0, 31.5, 40.0, 50.0, 63.0,
                                80.0, 100.0, 125.0, 160.0, 200.0, 250.0, 315.0,
                                400.0,500.0, 630.0, 800.0, 1000.0, 1250.0, 1600.0,
                                2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0,
                                8000.0, 10000.0, 12500.0, 16000.0, 20000.0])

    centerFreqs_oct = np.array([16.0, 31.5, 63.0, 125.0, 250.0, 500.0,
                                1000.0, 2000.0, 4000.0, 8000.0, 16000.0])

    filter_order = 2

    def __init__(self, data, sr):
        self.data = data
        self.output = []
        self.sr = sr
        self.centerFreqs = []

        return

    def band_filter(self, f_min, f_max, bandWidth = "third"):
        if f_min > f_max: raise RuntimeError("Minimum frequency greater than maximum frequency")
        if bandWidth == "third":
            centerFreqs = self.centerFreqs_3rd
            band_factor = 3
        elif bandWidth == "octave":
            centerFreqs = self.centerFreqs_oct
            band_factor = 1
        else:
            raise RuntimeError("Unsupported frequency bandWidth")
        f_min = min(centerFreqs, key=lambda x:abs(x-f_min))
        f_max = min(centerFreqs, key=lambda x:abs(x-f_max))
        f_min_idx = min(np.where(centerFreqs == f_min)).item()
        f_max_idx = min(np.where(centerFreqs == f_max)).item()
        centerFreqs = centerFreqs[f_min_idx:f_max_idx+1]
        for f_c in centerFreqs:
            f_low = f_c * 2 ** (-1 / ( 2 * band_factor ))
            f_high = f_c * 2 ** (1 / ( 2 * band_factor ))
            filtered_data = self.butter_bandpass_filter(f_low, f_high)
            self.output.append(filtered_data)
            self.centerFreqs = centerFreqs

        return self.output, self.centerFreqs

    def change_order(self, order):
        self.filter_order = order

        return

    def butter_bandpass(self, f_low, f_high, order=2):
        '''Creates the coeficients for a butterworth bandpass filter'''
        nyq = 0.5 * self.sr
        low = f_low / nyq
        high = f_high / nyq
        if order == 2:
            b, a = butter(self.filter_order, [low, high], btype='band')

            return b, a
        else:
            sos = butter(self.filter_order, [low, high], btype='band', output='sos')

            return sos


    def butter_bandpass_filter(self, f_low, f_high, order=2):
        '''Filters a signal using the filter created by butter_bandpass method'''
        if order == 2:
            b, a = self.butter_bandpass(f_low, f_high, order)
            # y = lfilter(b, a, self.data)
            y = filtfilt(b, a, self.data)

            return y
        else:
            sos = self.butter_bandpass(f_low, f_high, order)
            # y = sosfilt(sos, self.data)
            y = sosfiltfilt(sos, self.data)

            return y


    def removeDelay(self):
        IR_a = gz.a_weighting(self.data, self.sr)
        xcross_idx = np.where(np.diff(np.sign(IR_a)))[0]
        avg = np.empty(len(xcross_idx)-1)

        for i, _ in enumerate(xcross_idx[1:]):
            avg[i] = np.average(abs(self.data[xcross_idx[i]:xcross_idx[i+1]+1]))

        start_idx = xcross_idx[np.argmax(avg > 2e-4)]

        IR_head = self.data[:start_idx]
        IR_tail = self.data[start_idx:]

        return IR_head, IR_tail


    def find_nearest(array,value):
        idx = (np.abs(array-value)).argsort()[:2]

        return idx


    def fit_lorentz(p, x, y):
        '''Finds the coeficients for the Lorentz distribution that
        fits best to the data x,y'''
        return curve_fit(self.lorentz, x, y, p0 = p)


    def dB_clipp(self,threshold):
        '''Cuts a signal below a certain threshold'''
        Level = 20 * np.log10(abs(self.data))
        cut = np.where(Level > thershold)
        out = [Level[cut],self.data[cut],cut]

        return out


def findTau(data, tVec, threshold):
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
    testPoints = range(2000, len(data), 10)
    for testPoint in testPoints:
        b = (threshold - data[maxPoint])/(tVec[testPoint]-tVec[maxPoint])
        y2 = b*tVec[maxPoint:testPoint] + data[maxPoint]
        y = data[maxPoint:testPoint]

        rmsdiff.append(np.sqrt(((y - y2) ** 2).sum() / len(y)))
    tau = testPoints[np.argmin(rmsdiff)]

    return tau

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

def amp2dB(data):
    '''Converts amplitude to dB, using the maximum
    value as reference'''
    return 20 * np.log10((data)/max(data))

def lorentz(x, *p):
    '''Creates a Lorentz destribution using I, gamma and x0 values'''
    I, gamma, x0 = p

    return I * gamma**2 / ((x - x0)**2 + gamma**2)

