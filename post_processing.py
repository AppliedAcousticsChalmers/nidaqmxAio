#Python libraries
from plotly import tools
import plotly.offline as py
import plotly.graph_objs as go
# import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import signal, array
from tqdm import tqdm
#Program libraries
import utils as gz


def mic_calibration(measurements, sensitivity):
    #Reading the configuration
    keys = list(measurements.keys())
    sensitivity = float(sensitivity) * float(1e-3)
    calMeasurement = min([i1 for i1, s1 in enumerate(keys) if "cDAQ" in s1])
    values = list(measurements.values())
    data = np.array(values[calMeasurement])
    Lpmeas = max(gz.amp2db((data ) / (sensitivity) ,ref=2e-5))
    calibration = [10 **((93.9794 - Lpmeas)/20), sensitivity]

    return calibration

def h1_estimator_live(current_buffer, previous_buffer, blockidx, calibrationData, micAmp, refCH):
    current_buffer[:refCH,:] = (current_buffer[:refCH,:] / micAmp) * calibrationData[0]
    current_buffer[refCH+1:,:] = (current_buffer[refCH+1:,:] / micAmp) * calibrationData[0]
    blockSize = int(len(current_buffer[0]))
    nCHin = int(len(current_buffer))
    #window creation
    win = signal.hann(blockSize, sym=True)
    scaling = math.sqrt(np.sum(win ** 2) / blockSize)
    win =  win / scaling


    spectra = np.zeros((nCHin, blockSize), dtype=complex)
    for chidx in range(nCHin):
        #Applying Hanning window
        current_buffer[chidx,...] = current_buffer[chidx,...]*win
        #fft
        spectra[chidx,...] = np.fft.fft(current_buffer[chidx,...]) / blockSize
    #Calculating Sij, ASij, ADij (ij corresponds to each channel)
    S  = np.zeros((nCHin,nCHin,blockSize), dtype=complex)
    AD = np.zeros((nCHin,nCHin,blockSize), dtype=complex)
    AS = np.zeros((nCHin,nCHin,int(blockSize // 2)), dtype=complex)
    H = np.zeros((nCHin,int(blockSize // 2)), dtype=complex)
    H_phase = np.zeros((nCHin,int(blockSize // 2)), dtype=float)
    HD = np.zeros((nCHin,int(blockSize)), dtype=complex)
    IR = np.zeros((nCHin,int(blockSize)), dtype=complex)
    HdB  = np.zeros((nCHin, int(blockSize // 2)))
    gamma2 = np.zeros((nCHin, int(blockSize // 2)))
    for i in range(nCHin):
        for j in range(nCHin):
            if i!=j and i!=refCH: continue
            S[i,j,...] = np.multiply(spectra.conj()[i,...], spectra[j,...])
            if blockidx == 0:
                AD[i,j,...] = S[i,j,...]
            else:
                AD[i,j,...] = previous_buffer[i,j,...] - ((previous_buffer[i,j,...] - S[i,j,...]) / (blockidx + 1))
            AS[i,j,0] = AD[i,j,0]
            AS[i,j,-1] = AD[i,j,int(blockSize) // 2]
            AS[i,j,1:-1] = 2 * AD[i,j,1:int(blockSize // 2) - 1]
    for i in range(nCHin):
        if  i==refCH: continue
        H[i,...] = np.divide(AS[refCH,i,...], AS.real[refCH,refCH,...])
        HD[i,...] = np.divide(AD[refCH,i,...], AD[refCH,refCH,...])
        IR[i,...] = np.fft.ifft(HD[i,...])
        H_phase[i,...] = np.unwrap(np.angle(H[i,...])) #* np.pi/ sr
        gamma2[i,...] = (abs(AS[refCH,i,...] ** 2)) / (np.multiply(AS.real[refCH,refCH,...], AS.real[i,i,...]))
        HdB[i,...] = gz.amp2db(H[i,...])
    previous_buffer = AD
    blockidx += 1
    return previous_buffer, spectra, blockidx, H, HdB, H_phase, IR, gamma2

def h1_estimator(filename, blockSize, calibrationData):

    #Reading the configuration
    print("Importing raw data... ", end="")
    inputfile = np.load(filename)
    measurements = inputfile.item()
    measuredFRFs = measurements
    simTime = measurements["simulationTime"]
    sr = measurements["SampleRate"]
    micAmp = measurements["micAmp"]
    blockSize = int(blockSize) #is gonna be handled as a command line input
    signalUnpadded = measurements["Unpadded_signal"]
    #Initial calculations
    keys = list(measurements.keys())
    dataCH = [i1 for i1, s1 in enumerate(keys) if "cDAQ" in s1]
    firstCH = dataCH[0]
    lastCH = dataCH[-1]
    keys = keys[firstCH:lastCH+1]
    refCH = [i2 for i2, s2 in enumerate(keys) if measurements['Reference_channel'] in s2][0]
    print(refCH)
    values = list(measurements.values())
    data = np.array(values[firstCH:lastCH+1])
    nCHin = int(len(data))
    print("Completed")
    # Removing the zero padding of the signal
    if simTime < 5:
        print("Removing the zero padding of the signal (This will take a while)... ", end="")
        convxymin = np.convolve((data[int(refCH), ...]), np.flipud(signalUnpadded), 'valid')
        convxymax = np.convolve(np.flipud(data[int(refCH), ...]), signalUnpadded, 'valid')

        pad_start = np.argmax(convxymin)
        pad_end = np.argmax(convxymax)
        data = (data[..., pad_start:-pad_end] / calibrationData[1] * calibrationData[0]) / micAmp
        print("Completed")
    else:
        data[:refCH,:] = (data[:refCH,:] / micAmp) * calibrationData[0]
        data[refCH+1:,:] = (data[refCH+1:,:] / micAmp) * calibrationData[0]

    numberOfsamples = int(len(data[0]))
    #Spectra calculations
    print("Calculating the spectra... ", end="")
    zeropad = blockSize - numberOfsamples % blockSize
    nBlocks = int(numberOfsamples // blockSize)
    data = np.pad(data, [[0, 0], [0, int(zeropad)]], 'constant', constant_values = 0)
    dataMatrix = np.zeros((nCHin, blockSize, nBlocks))
    spectraMatrix = np.zeros((nCHin, blockSize, nBlocks), dtype=complex)
    print("Completed")


    #window creation
    print("Applying Hanning window... ", end="")
    win = signal.hann(blockSize, sym=True)
    scaling = math.sqrt(np.sum(win ** 2) / blockSize)
    win = win / scaling
    print("Completed")


    #Arranging the data into blocks, Applying window and calculating the fft
    print("Arranging the windowed spectra into blocks")
    for chidx in range(nCHin):
        i0 = 0
        pbar1 = tqdm(total=nBlocks)
        for idx in range(nBlocks):
            #Arranging into blocks, applying the window
            dataMatrix[chidx, :, idx] = data[chidx, i0:i0+blockSize] * win
            #Calculating the fft
            spectraMatrix[chidx, :, idx]= np.fft.fft(dataMatrix[chidx, :, idx]) / blockSize
            i0 += blockSize
            pbar1.update()
        pbar1.close()


    #Calculating Sij, ASij, AGij, Hij, gamma2ij (ij corresponds to each channel)
    print("Calculating the Double sided auto and cross spectra")
    S = np.zeros((nCHin,nCHin,blockSize), dtype=complex)
    AD = np.zeros((nCHin,nCHin,blockSize), dtype=complex)
    AS = np.zeros((nCHin,nCHin,blockSize // 2), dtype=complex)

    for i in range(nCHin):
        for j in range(nCHin):
            if i!=j and i!=refCH: continue
            pbar2 = tqdm(total=nBlocks)
            for blockidx in range(nBlocks):
                S[i,j,:] = np.multiply(spectraMatrix.conj()[i,:,blockidx], spectraMatrix[j,:,blockidx])
                if blockidx == 0:
                    AD[i,j,:] = S[i,j,:]
                else:
                    AD[i,j,:] -= (AD[i,j,:] - S[i,j,:]) / (blockidx + 1)
                pbar2.update()
            AS[i,j,0] = AD[i,j,0]
            AS[i,j,-1] = AD[i,j,int(blockSize // 2)]
            AS[i,j,1:-1] = 2 * AD[i,j,1:int(blockSize // 2 - 1)]
            pbar2.close()

    #Calculating ASij, Hij, gamma2ij
    print("Calculating the single sided spectra, FRFs and coherence")
    H = np.zeros((nCHin,blockSize // 2), dtype=complex)
    H_phase = np.zeros((nCHin,blockSize // 2), dtype=float)
    HD = np.zeros((nCHin,blockSize), dtype=complex)
    IR = np.zeros((nCHin,blockSize), dtype=complex)
    HdB = np.zeros((nCHin,blockSize // 2), dtype=float)
    gamma2 = np.zeros((nCHin,blockSize // 2), dtype=float)
    SNR = np.zeros((nCHin,blockSize // 2), dtype=float)

    pbar3 = tqdm(total=nCHin)
    for i in range(nCHin):
        if i==refCH: pbar3.update();continue
        H[i,:] = np.divide(AS[refCH,i,:], AS.real[refCH,refCH,:])
        HD[i,:] = np.divide(AD[refCH,i,:], AD.real[refCH,refCH,:])
        IR[i,:] = np.fft.ifft(HD[i,:])
        H_phase[i,:] = np.unwrap(np.angle(H[i,:])) #* np.pi/ sr
        gamma2[i,:] = (abs(AS[refCH,i,:] ** 2)) / (np.multiply(AS.real[refCH,refCH,:], AS.real[i,i,:]))
        SNR[i,:] = gamma2[i,:] / (1 - gamma2[i,:])
        HdB[i,:] = gz.amp2db(H[i,:], ref=1) #- 20 * np.log10(abs(micAmp))
        pbar3.update()
    pbar3.close()

    print("Saving the processed data... ", end="")
    fftfreq = np.fft.rfftfreq(blockSize, 1 / sr)
    tVec = np.linspace(0, len(IR[0,:]) / sr, len(IR[0,:]))
    measuredFRFs.update({"Hij":H, "HdBij":HdB, "HDij":HD, "IRij":IR, "gamma2ij":gamma2, "fftfreq":fftfreq, "tVec":tVec})
    print("Completed")

    for pltidx in range(nCHin):
        if pltidx == refCH: continue
        plot1y = HdB[pltidx,:]
        plot2y = gamma2[pltidx,:]
        plot3y = IR.real[pltidx,:]
        plot4y = H_phase[pltidx,:]
        freqs = fftfreq[:-1]
        tvec = tVec

    #     ax1.set_ylabel(r"$H1_{[%i,%i]}$ in dB"%(refCH[0],pltidx))
        trace1 = go.Scatter(
            x=freqs,
            y=plot1y,
            name='$H_{[%i,%i]}$ in dB'%(refCH,pltidx)
        )

        trace2 = go.Scatter(
            x=freqs,
            y=plot2y,
            name='$\\gamma^2_{[%i,%i]}$'%(refCH,pltidx)
        )

        trace3 = go.Scatter(
            x=tvec,
            y=plot3y,
            name='$\\text{IR}_{[%i,%i]}$'%(refCH,pltidx)
        )

        trace4 = go.Scatter(
            x=freqs,
            y=plot4y,
            name='$Phase_{H,[%i,%i]}$'%(refCH,pltidx)
        )

        fig = tools.make_subplots(rows=2, cols=2, subplot_titles=('HdB', 'gamma2', 'IR', 'H_phase'))

        fig.append_trace(trace1, 1, 1)
        fig.append_trace(trace2, 1, 2)
        fig.append_trace(trace3, 2, 1)
        fig.append_trace(trace4, 2, 2)

        fig['layout'].update(title=filename,
                             xaxis1=dict(title='frequency in Hz', type='log', autorange=True),
                             xaxis2=dict(title='frequency in Hz', type='log', autorange=True),
                             xaxis3=dict(title='time in s', type='lin', autorange=True),
                             xaxis4=dict(title='frequency in Hz', type='log', autorange=True)
        )
        py.plot(fig, filename=filename +'.html')
    print("Process Finished")


    return measuredFRFs

def calT60(filename):

    #Reading the configuration
    print("Importing IR data... ", end="")
    inputfile = np.load(filename)
    measurements = inputfile.item()
    measuredFRFs = measurements
    # simTime = measurements["simulationTime"]
    sr = measurements["SampleRate"]
    #Initial calculations
    keys = list(measurements.keys())
    dataCH = [i1 for i1, s1 in enumerate(keys) if "cDAQ" in s1]
    firstCH = dataCH[0]
    lastCH = dataCH[-1]
    keys = keys[firstCH:lastCH+1]
    refCH = [i2 for i2, s2 in enumerate(keys) if measurements['Reference_channel'] in s2][0]
    values = list(measurements.values())
    print(values, refCH)
    del values[refCH]
    print(values)
    measuredFRFs.update({"Input_Channel_names":keys})
    data = np.array(values[firstCH:lastCH+1])
    nCHin = int(len(data))
    print("Completed")
