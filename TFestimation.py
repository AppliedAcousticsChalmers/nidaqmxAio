# Python libraries
import numpy as np
from scipy import signal, array
import math

# Program libraries
import filtersAndMathtools as flm

def h1_estimator(current_buffer, previous_buffer, blockidx, calibrationData, micAmp, refCH):
    '''
    Calculates the TF of a system using the H1 estimator method

    The function works in a block by block fashion, averaging the current input with
    it's previous output while keeping track the number of blocks that have been processed.
    For this reason it's usable also for real time processing of TF measurements.

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
    - H: the resulting single sided TF
    - HD: the resulting double sided TF
    - HdB: the single sided TF in dB with reference 1V
    - H_phase: the phase of the TF in dB reference 1V
    - IR: the resulting IR of the system, in V
    - gamma2: the resulting coherence
    '''

    #Applying calibration data
    current_buffer[:refCH, :] = (current_buffer[:refCH, :] / micAmp) * calibrationData[0]
    current_buffer[refCH+1:, :] = (current_buffer[refCH+1:, :] / micAmp) * calibrationData[0]

    #Calculating the signal length and the number of channels
    blockSize = int(len(current_buffer[0]))
    rfftBlock = int(blockSize / 2 + 1)
    nCHin = int(len(current_buffer))

    # window creation
    win = signal.hann(blockSize, sym=True)
    scaling = math.sqrt(np.sum(win ** 2) / blockSize)
    win = win / scaling

    spectra = np.zeros((nCHin, blockSize), dtype=complex)
    for chidx in range(nCHin):
        # Applying Hanning window
        current_buffer[chidx, :] = current_buffer[chidx, :]*win
        # fft
        spectra[chidx, :] = np.fft.fft(current_buffer[chidx, :]) / blockSize

    # Calculating S, AS, AD
    S = np.zeros((nCHin, nCHin, blockSize), dtype=complex)
    AD = np.zeros((nCHin, nCHin, blockSize), dtype=complex)
    AS = np.zeros((nCHin, nCHin, rfftBlock), dtype=complex)
    H = np.zeros((nCHin, rfftBlock), dtype=complex)
    H_phase = np.zeros((nCHin, rfftBlock), dtype=float)
    HD = np.zeros((nCHin, int(blockSize)), dtype=complex)
    IR = np.zeros((nCHin, int(blockSize)), dtype=float)
    HdB = np.zeros((nCHin, rfftBlock), dtype=float)
    gamma2 = np.zeros((nCHin, rfftBlock), dtype=float)
    SNR = np.zeros((nCHin, rfftBlock), dtype=float)

    for i in range(nCHin):
        for j in range(nCHin):
            if i != j and i != refCH:
                continue
            S[i, j, :] = np.multiply(spectra.conj()[i, :], spectra[j, :])
            if blockidx == 0:
                AD[i, j, :] = S[i, j, :]
            else:
                AD[i, j, :] = previous_buffer[i, j, :] - ((previous_buffer[i, j, :] - S[i, j, :]) / (blockidx + 1))
            AS[i, j, 0] = AD[i, j, 0]
            AS[i, j, -1] = AD[i, j, rfftBlock]
            AS[i, j, 1:-1] = 2 * AD[i, j, 1:rfftBlock - 1]

    #Calculating H, HD, IR, H_phase, gamma2, HdB, SNR
    for i in range(nCHin):
        if i == refCH:
            continue
        H[i, :] = np.divide(AS[refCH, i, :], AS[refCH, refCH, :])
        HD[i, :] = np.divide(AD[refCH, i, :], AD[refCH, refCH, :])
        IR[i, :] = np.fft.ifft(HD[i, :]).real
        H_phase[i, :] = np.unwrap(np.angle(H[i, :]))
        gamma2[i, :] = (abs(AS[refCH, i, :] ** 2)) / (np.multiply(AS.real[refCH, refCH, :], AS.real[i, i, :]))
        HdB[i, :] = flm.amp2db(H[i, :])
        SNR[i, :] = gamma2[i, :] / (1 - gamma2[i, :])
    previous_buffer = AD
    blockidx += 1

    return previous_buffer, spectra, blockidx, H, HD, HdB, H_phase, IR, gamma2, SNR

def deconvolution(data, sr, blockSize, calibrationData, micAmp, refCH, f0, f1):
    '''
    Calculates the IR and TF using the deconvolution method.

    -INPUTS-
    - data: raw multichannel data in the form of a numpy array (number of channels, samples)
    - sr: sampling rate
    - blockSize: blockSize for the spectrogram's fft analysis and number of samples of the
    resulting IR that contain valuable information.
    - calibrationData: list with [ calibration coefficient, microphone sensitivity]
    - micAmp: microphone preamplificaltion
    - refCH: the index that corresponds to the reference channel in the data
    - f0: start frequency of the sweep
    - f1: stop frequency of the sweep

    -OUTPUTS-
    - H: the resulting single sided TF
    - HD: the resulting double sided TF
    - HdB: the single sided TF in dB with reference 1V
    - H_phase: the pahse of the TF in dB reference 1V
    - IR: the resulting IR of the system, in V
    - fftfreq: the corresponding frequency vector
    - tVec: the corresponding time vector
    - spectrogram: the resulting spectrogram. list:[time vector, frequency vector, magnitude]


    For optimal accuracy three conditions should be taken under consideration.

    1. The margins of the sweep should be larger that the required frequency range
    of measurement, in order to ensure that the filtering will be applied in a range
    where no information will be lost. Since this algorithm is created with measurements
    of room acoustics in mind, the applied filter is a 2nd order butterworth in order for
    it's IR to be short enough to be applied in a rooms IR, without losing information.

    2. Zero padding in the beginning of the signal should be used in order for the sweep
    to start from zero, avoiding non linearities in the fft analysis.

    3. Sufficient measurement time should be considered for better SNR, while the decay of
    the system should be also recorded after the excitation has stopped. This can be set
    using the cutoffTime parameter of the program.

    The calculation procedure is as follows.

    a. The calibrationData and the external microphone amplification is applied to the
    raw data signals

    b. The frequency vector of the fft analysis is calculated and the frequencies that
    contain meaningful information are detected. This refers to the frequencies of the
    fft analysis that are inside the frequency range of the sweep. Outside of this range
    the fft analysis will just contain numerical errors that are not useful, or posses
    any physical meaning. Therefore the spectra outside of the meaningful frequency range
    are set to zero. This procedure is done in the single sided spectrum, so two frequencies
    specify the useful frequency range, f0, f1, where f0 and f1 correspond to
    the beginning and the end of the sweep respectively.

    c. The spectra of the measurement channels are divided by the spectra of the reference
    channel, yielding the TF of the system. The initial IR is calculated using numpy's irfft.

    d. The initial IR is filtered using a 2nd order bandpass butterworth filter, with critical
    frequencies [ 2 * f0, f1 / 2].

    e. The resulting IR will be as long as the measurement, while distortions will be pushed toward
    it's end. Therefore it's TF will be quite noisy. For that reason the IR is cut to the part containing
    useful information, using the blockSize parameter.

    f. The TF and phase can be calculated using the appropriate fft analysis in the resulting IR.

    '''

    #Applying calibration Data

    data[:refCH, :] = (data[:refCH, :] / micAmp) * calibrationData[0]
    data[refCH+1:, :] = (data[refCH+1:, :] / micAmp) * calibrationData[0]

    #Calculating the signal length and the number of channels
    signalLength = int(len(data[0]))
    # rfft keeps only the real frequency bins of the fft plus the nyquist frequency and the DC
    rfftBlock_temp = int(signalLength / 2 + 1) #Size prior to IR end clipping
    rfftBlock = int(blockSize / 2 + 1) #Size after IR clipping
    nCHin = int(len(data))

    #Calculating the initial frequency vector
    fftfreq_temp = np.fft.rfftfreq(signalLength, 1 / sr)

    #Calculating the frequencies that have meaningful information
    first_rfreq = np.min(np.where( fftfreq_temp >= f0 ))
    last_rfreq = np.max(np.where( fftfreq_temp <= f1 ))

    spectrogramm = []
    spectra = np.zeros((nCHin, rfftBlock_temp), dtype=complex)

    for chidx in range(nCHin):
        # fft analysis
        spectra[chidx, ...] = np.fft.rfft(data[chidx, ...])

        # Spectrogram
        f_sp, t_sp, Sxx_sp = signal.spectrogram(data[chidx, ...], sr, nperseg=blockSize, noverlap= blockSize // 2, window='triang', scaling='spectrum', mode='magnitude' )
        Sxx_sp = 10 * np.log10(abs(Sxx_sp))
        spectrogramm.append([f_sp, t_sp, Sxx_sp])

    # Vector preallocation
    IR = np.zeros((nCHin, blockSize), dtype=float)
    HD = np.zeros((nCHin, rfftBlock), dtype=complex)
    H = np.zeros((nCHin,  rfftBlock), dtype=complex)
    HdB = np.zeros((nCHin, rfftBlock), dtype=float)
    H_phase = np.zeros((nCHin, rfftBlock), dtype=float)
    IR_temp = np.zeros((nCHin,  signalLength), dtype=float)
    HD_temp = np.zeros((nCHin, rfftBlock_temp), dtype=complex)

    # Loop over the channels
    for i in range(nCHin):
        if i == refCH:
            continue
        # Discarding frequencies outside the sweep range
        spectra[i, :first_rfreq] = 0
        spectra[i, last_rfreq:] = 0

        # Deconvolution
        HD_temp[i,:] = spectra[i,:] / spectra[refCH,:]
        IR_temp[i,:] = np.fft.irfft(HD_temp[i,:])

        # Filtering
        IR_temp[i,:]= flm.butterWorthFilter(IR_temp[i,:], sr, frequecies=[f0 * 2, f1 / 2], order=2, filter_type='bandpass', method='ba')

        # Cutting the the tail of the IR
        win = flm.half_hann(int(sr*0.01), 'right')
        IR[i,:] = IR_temp[i,:blockSize]
        IR[i,-len(win):] *= win

        # Final results
        HD[i,:] = np.fft.rfft(IR[i,:])
        H[i,:] = HD[i,:]
        HdB[i,:] = flm.amp2db(H[i,:])
        H_phase[i,:] = np.unwrap(np.angle(H[i,:]))

    tVec = np.linspace(0, blockSize / sr, blockSize)
    fftfreq = np.fft.rfftfreq(blockSize, 1 / sr)

    return H, HD, HdB, H_phase, IR, fftfreq, tVec, spectrogramm
