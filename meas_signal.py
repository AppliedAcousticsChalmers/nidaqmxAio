import numpy as np
import colorednoise as cn
from scipy import signal


def create_signal(sample_rate, time_of_signal, pad_samples, signal_type, voltage_out):
    """Create the output signal to be measured.

    :sample_rate: sample rate for signal generation (int)
    :time_of_signal: signal length in seconds (float/int)
    :padN: zero padding in front of the signal (int)
    :signal_type: curently "noise_pink", "noise_white", "tone", "sweep_logarithmic" etc.
    :voltage_out: peak output in voltage

    :returns: (signal, unpaddedSignal)

    """

    number_of_samples = int(time_of_signal * sample_rate)
    time_of_signal = number_of_samples / sample_rate
    time_vector = np.linspace(0, time_of_signal, number_of_samples)
    signal_unpadded = np.empty(len(time_vector))
    sig = signal_type[0].split("_")
    print(np.shape(signal_unpadded))
    if sig[0] == "noise":
        if sig[1] == "pink":
            signal_unpadded = cn.powerlaw_psd_gaussian(1, number_of_samples)
        elif sig[1] == "white":
            signal_unpadded = cn.powerlaw_psd_gaussian(0, number_of_samples)
        else:
            print('Non upported noise type. Reverting to pink noise')
            signal_unpadded = cn.powerlaw_psd_gaussian(1, number_of_samples)
    elif sig[0] == "tone":
        f = int(signal_type[1])
        signal_unpadded = np.sin(f * 2 * np.pi * time_vector)
    elif sig[0] == "sweep":
        method = sig[1]
        f0 = int(signal_type[1])
        f1 = int(signal_type[2])
        signal_unpadded = signal.chirp(time_vector, f0, time_of_signal, f1, method)
    else:
        print('Unlnown signal type. Reverting to pink noise')
        signal_unpadded = cn.powerlaw_psd_gaussian(1, number_of_samples)

    signal_unpadded /= np.max(abs(signal_unpadded)) / voltage_out
    signal_padded = np.pad(signal_unpadded, (pad_samples, 0), 'constant', constant_values=[0, 0])

    return signal_padded, signal_unpadded
