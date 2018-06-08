# Python libraries
import numpy as np
from scipy import signal, array
import math

import matplotlib.pyplot as plt
# Program libraries
import utils as gz
import filters as fl

def h1_estimator(current_buffer, previous_buffer, blockidx, calibrationData, micAmp, refCH):
    '''
    Calculates the TF of a system using the H1 estimator method

    The function works in a block by block fashion, averaging the current input with
    it's previous output while keepping track the number of blocks that have been processed.
    For this reason it's usable also for real time prossesing of TF measurements.

    -INPUTS-
    - current_buffer: current block of raw data in a numpy array (number of channels, samples)
    - previous_buffer: the previous output of the averaged double sided auto and cross spectra.
    A numpy array with shape ( number of channels, samples). Should be initialized with zeros
    for the first block of data.
    - blockidx: the number of blocks processed before the current. Updates internally.Should be initialized
    to zero.
    - calibrationData: list with [ calibration coefficient, microphone sensitivity]
    - micAmp: microphone preamplificaltion
    - refCH: the index that corresponds to the reference channel in the data

    -OUTPUTS-
    - previous_buffer: as in the inputs
    - spectra: the single sided spectra of the current buffer
    - blockidx: as in inputs
    - H: the resulting single sidded TF
    - HD: the resulting double sidded TF
    - HdB: the single sidded TF in dB with reference 1V
    - H_phase: the pahse of the TF in dB reference 1V
    - IR: the resulting IR of the system, in V
    - gamma2: the resulting coherence
    '''

    #Applying calibration data
    current_buffer[:refCH, :] = (current_buffer[:refCH, :] / micAmp) * calibrationData[0]
    current_buffer[refCH+1:, :] = (current_buffer[refCH+1:, :] / micAmp) * calibrationData[0]

    #Calculating the signal length and the number of channels
    blockSize = int(len(current_buffer[0]))
    nCHin = int(len(current_buffer))

    # window creation
    win = signal.hann(blockSize, sym=True)
    scaling = math.sqrt(np.sum(win ** 2) / blockSize)
    win = win / scaling

    spectra = np.zeros((nCHin, blockSize), dtype=complex)
    for chidx in range(nCHin):
        # Applying Hanning window
        current_buffer[chidx, ...] = current_buffer[chidx, ...]*win
        # fft
        spectra[chidx, ...] = np.fft.fft(current_buffer[chidx, ...]) / blockSize

    # Calculating S, AS, AD
    S = np.zeros((nCHin, nCHin, blockSize), dtype=complex)
    AD = np.zeros((nCHin, nCHin, blockSize), dtype=complex)
    AS = np.zeros((nCHin, nCHin, int(blockSize // 2)), dtype=complex)
    H = np.zeros((nCHin, int(blockSize // 2)), dtype=complex)
    H_phase = np.zeros((nCHin, int(blockSize // 2)), dtype=float)
    HD = np.zeros((nCHin, int(blockSize)), dtype=complex)
    IR = np.zeros((nCHin, int(blockSize)), dtype=complex)
    HdB = np.zeros((nCHin, int(blockSize // 2)))
    gamma2 = np.zeros((nCHin, int(blockSize // 2)))
    # SNR = np.zeros((nCHin, blockSize // 2), dtype=float)
    for i in range(nCHin):
        for j in range(nCHin):
            if i != j and i != refCH:
                continue
            S[i, j, ...] = np.multiply(spectra.conj()[i, ...], spectra[j, ...])
            if blockidx == 0:
                AD[i, j, ...] = S[i, j, ...]
            else:
                AD[i, j, ...] = previous_buffer[i, j, ...] - ((previous_buffer[i, j, ...] - S[i, j, ...]) / (blockidx + 1))
            AS[i, j, 0] = AD[i, j, 0]
            AS[i, j, -1] = AD[i, j, int(blockSize) // 2]
            AS[i, j, 1:-1] = 2 * AD[i, j, 1:int(blockSize // 2) - 1]

    #Calculating H, HD, IR, H_phase, gamma2, HdB, SNR
    for i in range(nCHin):
        if i == refCH:
            continue
        H[i, ...] = np.divide(AS[refCH, i, ...], AS.real[refCH, refCH, ...])
        HD[i, ...] = np.divide(AD[refCH, i, ...], AD[refCH, refCH, ...])
        IR[i, ...] = np.fft.ifft(HD[i, ...])
        H_phase[i, ...] = np.unwrap(np.angle(H[i, ...]))
        gamma2[i, ...] = (abs(AS[refCH, i, ...] ** 2)) / (np.multiply(AS.real[refCH, refCH, ...], AS.real[i, i, ...]))
        HdB[i, ...] = gz.amp2db(H[i, ...])
        # SNR[i, :] = gamma2[i, :] / (1 - gamma2[i, :])
    previous_buffer = AD
    blockidx += 1

    return previous_buffer, spectra, blockidx, H, HD, HdB, H_phase, IR, gamma2 #,SNR



def deconvolution(data, sr, blockSize, calibrationData, micAmp, refCH, f0, f1):
    '''
    Calculates the IR and TF using the deconvolution method.

    -INPUTS-
    - data: raw multychannel data in the form of a numpy array (number of channels, samples)
    - sr: sampling rate
    - blockSize: blockSize for the spectrogram's fft analysis and number of samples of the
    resulting IR that contain valuable information.
    - calibrationData: list with [ calibration coefficient, microphone sensitivity]
    - micAmp: microphone preamplificaltion
    - refCH: the index that corresponds to the reference channel in the data
    - f0: start frequecy of the sweep
    - f1: stop frequency of the sweep

    -OUTPUTS-
    - H: the resulting single sidded TF
    - HD: the resulting double sidded TF
    - HdB: the single sidded TF in dB with reference 1V
    - H_phase: the pahse of the TF in dB reference 1V
    - IR: the resulting IR of the system, in V
    - fftfreq: the corresponding frequency vector
    - tVec: the corresponding time vector
    - spectrogramm: the resulting spectrogramm. list:[time vector, frequency vector, magnitude]


    For optimal accurasy three conditions should be taken under consideration.

    1. The margins of the sweep should be larger that the required frequency range
    of measurement, in order to ensure that the filtering will be applied in a range
    where no information will be lost. Since this algorithm is created with measurements
    of room acoustics in mind, the applied filter is a 1st order butterworth in order for
    it's IR to be short enough to be applied in a rooms IR, without losing information.

    2. Zeropadding in the beginning of the signal should be used in order for the sweep
    to start from zero, avoiding non linearities in the fft analysis.

    3. Sufficient measurement time should be considered for better SNR, while the decay of
    the system should be also recorded after the excitation has stopped. This can be set
    using the cutoffTime parameter of the program.

    The calculation procedure is as follows.

    a. The calibrationData and the external microphone amplification is applied to the
    raw data signals

    b. The frequency vector of the fft analisys is calulated and the frequencies that
    contain meaningfull information are detected. This refers to the frequencies of the
    fft analisys that are inside the frequency range of the sweep. Outside of this range
    the fft analisys will just contain numerical errors that are not usefull, or posses
    any physical meaning. Therefore the spectra outside of the meanigfull frequency range
    are set to zero. This procedure is done in the double sided spectrum, so four frequncies
    specyfy the usefull frequency range, f0, f1, -f1, -f1, where f0 and f1 correspond to
    the beginning and the end of the sweep respectively.

    c. The spectra of the measurement channels are divided by the spectra of the the reference
    channel, yielding the TF of the system. The initial IR is calculated using the ifft of the TFcalc

    d. The initial IR is filtered using a 1st order bandpass butterworth filter, with critical
    frequencies [ 2 * f0, f1 / 2].

    e. The resulting IR will be as long as the measurement and the therefore it's TF will be quite
    noisy. For that reason the IR is cutt to the part which it contains usefull information, using
    the blockSize parameter.

    f. The TF and phase can be calculated using the apprpriate ft analisys in the resulting IR.

    '''

    #Applying calibration Data

    data[:refCH, :] = (data[:refCH, :] / micAmp) * calibrationData[0]
    data[refCH+1:, :] = (data[refCH+1:, :] / micAmp) * calibrationData[0]

    #Calculating the signal length and the number of channels
    signalLength = int(len(data[0]))
    nCHin = int(len(data))

    #Calculating the initial frequency vector
    fftfreq = np.fft.fftfreq(signalLength, 1 / sr)

    #Calculating the frequencies that have meaningfull information
    first_rfreq = np.min(np.where( fftfreq[:len(fftfreq)//2] >= f0 ))
    last_rfreq = np.max(np.where( fftfreq[:len(fftfreq)//2] <= f1 ))
    first_ifreq = len(fftfreq)//2 + np.min(np.where( fftfreq[len(fftfreq)//2:] >= -f0 ))
    last_ifreq =  len(fftfreq)//2 + np.max(np.where( fftfreq[len(fftfreq)//2:] <= -f1 ))

    spectrogramm = []
    spectra = np.zeros((nCHin, signalLength), dtype=complex)

    for chidx in range(nCHin):
        # fft analisys
        spectra[chidx, ...] = np.fft.fft(data[chidx, ...])

        # Spectrogram
        f_sp, t_sp, Sxx_sp = signal.spectrogram(data[chidx, ...], sr, nperseg=blockSize, noverlap= blockSize // 2, window='triang', scaling='density' )
        Sxx_sp = 20 *np.log10(abs(Sxx_sp))
        spectrogramm.append([f_sp, t_sp, Sxx_sp])

    # Vector preallocation
    HD = np.zeros((nCHin, blockSize), dtype=complex)
    H = np.zeros((nCHin,  blockSize//2), dtype=complex)
    HdB = np.zeros((nCHin, blockSize//2), dtype=float)
    H_phase = np.zeros((nCHin, blockSize//2), dtype=float)
    IR = np.zeros((nCHin, blockSize), dtype=complex)
    IR_temp = np.zeros((nCHin,  len(spectra[0])), dtype=complex)
    HD_temp = np.zeros((nCHin, len(spectra[0])), dtype=complex)

    # Loop over the channels
    for i in range(nCHin):
        if i == refCH:
            continue
        # Discarding frequencies outsided the sweep range
        spectra[i, :first_rfreq] = 0
        spectra[i, last_rfreq:first_ifreq] = 0
        spectra[i, last_ifreq:] = 0

        # Deconvolution
        HD_temp[i,:] = spectra[i,:] / spectra[refCH,:]
        IR_temp[i,:] = (np.fft.ifft(HD_temp[i,:]))

        # Filtering
        IR_temp[i,:] = fl.filters(IR_temp[i,:], sr).butter_bandpass_filter(f0 * 2 , f1 / 2 , 1)

        # Cutting the the tail of the IR
        IR[i,:] = IR_temp[i,:blockSize]

        # Final results
        HD[i,:] = np.fft.fft(IR[i,:])
        H[i,:] = HD[i,:len(HD[i,:])//2]
        HdB[i,:] = gz.amp2db(H[i,:])
        H_phase[i,:] = np.unwrap(np.angle(H[i,:]))

    tVec = np.linspace(0, len(IR[0, :]) / sr, len(IR[0, :]))
    fftfreq = np.fft.fftfreq(blockSize, 1 / sr)

    return H, HD, HdB, H_phase, IR, fftfreq, tVec, spectrogramm
