import numpy as np
import colorednoise as cn
from scipy import signal
import matplotlib.pyplot as plt
import scipy.io as io

#Program libraries
import filtersAndMathtools as flm


def create_signal(sample_rate, time_of_signal, pad_samples, signal_type, voltage_out, cutoffTime=0):
    """Create the output signal to be measured.

    :sample_rate: sample rate for signal generation (int)
    :time_of_signal: signal length in seconds (float/int)
    :padN: zero padding in front of the signal (int)
    :signal_type: curently "noise_pink", "noise_white", "tone", "sweep_logarithmic" etc.
    :voltage_out: peak output in voltage
    :cutoffTime: time in seconds that the signal is padded at the end. Used for interrupted noise measurements

    :returns: (signal, unpadded_Signal)

    """

    number_of_samples = int(time_of_signal * sample_rate)
    time_of_signal = number_of_samples / sample_rate
    time_vector = np.linspace(0, time_of_signal, number_of_samples)
    signal_unpadded = np.empty(len(time_vector))
    sig = signal_type[0].split("_")
    pad_end = int(cutoffTime * sample_rate)
    cutoffTime = pad_end / sample_rate
    if sig[0] == "noise":
        if sig[1] == "pink":
            signal_unpadded = cn.powerlaw_psd_gaussian(1, number_of_samples)
        elif sig[1] == "white":
            signal_unpadded = cn.powerlaw_psd_gaussian(0, number_of_samples)
        else:
            print('Non supported noise type. Reverting to pink noise')
            signal_unpadded = cn.powerlaw_psd_gaussian(1, number_of_samples)
        win_left = flm.half_hann(int(sample_rate * 0.001), 'left')
        win_right = flm.half_hann(int(sample_rate * 0.001), 'right')
        signal_unpadded[-len(win_right):] *= win_right
        signal_unpadded[:len(win_right)] *= win_left
    elif sig[0] == "tone":
        f = float(signal_type[1])
        win_right = flm.half_hann(int(sample_rate * 0.001), 'right')
        win_left = flm.half_hann(int(sample_rate * 0.001), 'left')
        signal_unpadded = np.sin(float(f) * 2 * np.pi * time_vector)
        signal_unpadded[-len(win_right):] *= win_right
        signal_unpadded[:len(win_right)] *= win_left
    elif sig[0] == "sweep":
        method = sig[1]
        f0 = int(signal_type[1])
        f1 = int(signal_type[2])
        win_right = flm.half_hann(int(sample_rate * 0.001), 'right')
        win_left = flm.half_hann(int(sample_rate * 0.001), 'left')
        signal_unpadded = signal.chirp(time_vector, f0, time_of_signal, f1, method, phi=270)
        signal_unpadded[-len(win_right):] *= win_right
        signal_unpadded[:len(win_left)] *= win_left
    elif sig[0] == "matLab":
        matFile = io.loadmat(signal_type[1]+'.mat')
        signal_unpadded = matFile['audio'][0]
        time_of_signal = int(len(signal_unpadded) / sample_rate)
    else:
        print('Unlnown signal type. Reverting to pink noise')
        signal_unpadded = cn.powerlaw_psd_gaussian(1, number_of_samples)

    signal_unpadded /= np.max(abs(signal_unpadded)) / voltage_out
    signal_padded = np.pad(signal_unpadded, (pad_samples, pad_end), 'constant', constant_values=[0, 0])
    signal_size_in_samples = len(signal_padded)

    return [signal_padded, signal_unpadded, signal_size_in_samples, time_of_signal, cutoffTime]

def testSig(sr, t, tones=[[100, 1, 0, 'sin']], noiseType='no', plotting=False):
    '''
    Creates signals composed of tones and noise for test purposes. This function is not used anywhere in the
    program.

    sr - sample rate
    t - time in seconds
    tones - list of tones that the signal will contain
    format: tones=[tone1,tone2,...], where tone1 = [frequency, Amplitute, starting phase, type(sin or cos)]
    noise - noise added to the resulting signal
    format: noise=[type, amplitude in percetage of the overall max amplitude of the tonal part of the signal]
    plotting - plots the resulting signal

    '''
    tVec = np.linspace(0, t, t * sr)
    number_of_samples = len(tVec)
    testSig = 0
    for i in range(0, len(tones)):
        A = tones[i][1]
        f = tones[i][0]
        if len(tones[i]) == 2 :
            phi = 0
            tone_type = 'sin'
        elif len(tones[i]) == 3 :
            tone_type = 'sin'
        else:
            phi = tones[i][2]
            tone_type = tones[i][3]
        if tone_type == 'sin':
            testSig += A * np.sin(2 * np.pi * f * tVec + phi)
        elif tone_type == 'cos':
            testSig += A * np.cos(2 * np.pi * f * tVec + phi)
        else:
            print("Unknown tone type skipping...")
    #Adding noise
    if noiseType == 'no':
        noise = 0
    else:
        if noiseType[0] == 'pink':
            noise = cn.powerlaw_psd_gaussian(1, number_of_samples)
        elif noiseType[0] == 'white':
            noise = cn.powerlaw_psd_gaussian(0, number_of_samples)
        else:
            noise = np.random.randn(number_of_samples)
        noise /= np.max(abs(testSig))
        noise *= noiseType[1]
    testSig += noise

    if plotting == True:
        plt.plot(tVec, testSig)
        plt.show()

    return testSig, tVec, tones, noise
