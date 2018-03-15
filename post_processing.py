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
    calMeasurement = min([i1 for i1, s1 in enumerate(keys) if "cDAQ" in s1])
    values = list(measurements.values())
    data = np.array(values[calMeasurement])
    Lpmeas = max(gz.amp2db((data ) / (float(sensitivity) * float(1e-3)),ref=2e-5))
    calibration = [10 **((93.9794 - Lpmeas)/20), sensitivity]

    return calibration

def h1_estimator_live(current_buffer, previous_buffer, blockidx, calibrationData, micAmp, refCH):
    current_buffer = ((np.array(current_buffer) / calibrationData[1]) * calibrationData[0]) / micAmp
    # previous_buffer = np.array(previous_buffer)
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
            # S[i,j,0:2] = S[i,j,0:2]*0 + S[i,j,3]
            # S[i,j,-2:] = S[i,j,0:1]*0 + S[i,j,3]
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
    # keys = keys[dataCH[0]:dataCH[-1]+1]
    firstCH = dataCH[0]
    lastCH = dataCH[-1]
    keys = keys[firstCH:lastCH+1]
    refCH = [i2 for i2, s2 in enumerate(keys) if measurements['Reference_channel'] in s2]
    values = list(measurements.values())
    measuredFRFs.update({"Input_Channel_names":keys})
    data = np.array(values[firstCH:lastCH+1])
    nCHin = int(len(data))
    print("Completed")

    # Removing the zero padding of the signal
    if simTime < 5:
        print("Removing the zero padding of the signal (This will take a while)... ", end="")
        convxymin = np.convolve((data[int(refCH[0]), ...]), np.flipud(signalUnpadded), 'valid')
        convxymax = np.convolve(np.flipud(data[int(refCH[0]), ...]), signalUnpadded, 'valid')

        pad_start = np.argmax(convxymin)
        pad_end = np.argmax(convxymax)
        data = (data[..., pad_start:-pad_end] / calibrationData[1] * calibrationData[0]) / micAmp
        print("Completed")
    else:
        data = (data / calibrationData[1] * calibrationData[0]) / micAmp

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
            if i!=j and i!=refCH[0]: continue
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
    H = np.zeros((nCHin,nCHin,blockSize // 2), dtype=complex)
    H_phase = np.zeros((nCHin,nCHin,blockSize // 2), dtype=float)
    HD = np.zeros((nCHin,nCHin,blockSize), dtype=complex)
    IR = np.zeros((nCHin,nCHin,blockSize), dtype=complex)
    HdB = np.zeros((nCHin,nCHin,blockSize // 2), dtype=float)
    gamma2 = np.zeros((nCHin,nCHin,blockSize // 2), dtype=float)
    SNR = np.zeros((nCHin,nCHin,blockSize // 2), dtype=float)

    pbar3 = tqdm(total=nCHin)
    for i in range(nCHin):
        for j in range(nCHin):
            if i!=j and i!=refCH[0]: continue
            H[i,j,:] = np.divide(AS[i,j,:], AS.real[i,i,:])
            HD[i,j,:] = np.divide(AD[i,j,:], AD.real[i,i,:])
            IR[i,j,:] = np.fft.ifft(HD[i,j,:])
            H_phase[i,j,:] = np.unwrap(np.angle(H[i,j,:])) #* np.pi/ sr
            gamma2[i,j,:] = (abs(AS[i,j,:] ** 2)) / (np.multiply(AS.real[i,i,:], AS.real[j,j,:]))
            SNR[i,j,:] = gamma2[i,j,:] / (1 - gamma2[i,j,:])
            HdB[i,j,:] = gz.amp2db(H[i,j,:], ref=2e-5 * 2 ** 0.5)
        pbar3.update()
    pbar3.close()

    print("Saving the processed data... ", end="")
    fftfreq = np.fft.rfftfreq(blockSize, 1 / sr)
    tVec = np.linspace(0, len(IR[0,0,:]) / sr, len(IR[0,0,:]))
    measuredFRFs.update({"Hij":H, "HdBij":HdB, "HDij":HD, "gamma2ij":gamma2, "fftfreq":fftfreq, "tVec":tVec})
    print("Completed")

    #Plotting
    # for pltidx in range(len(H)):
    #     if pltidx == refCH[0]: continue
    #     fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)

    #     ax1.semilogx(fftfreq[:-1], HdB[refCH[0],pltidx,...])
    #     ax1.set_xlabel("frequency in Hz")

    #     ax2.semilogx(fftfreq[:-1], gamma2[refCH[0],pltidx,...]) 
    #     # ax2.semilogx(fftfreq[:-1], SNR[refCH[0],pltidx,...])
    #     ax2.set_ylabel(r"$\gamma^2_{[%i,%i]},$"%(refCH[0],pltidx))
    #     ax2.set_xlabel("frequency in Hz")

    #     ax3.plot(tVec,IR.real[refCH[0],pltidx,...])
    #     ax3.set_ylabel(r"$IR_{[%i,%i]}$"%(refCH[0],pltidx))
    #     ax3.set_xlabel("time is s")

    #     ax4.plot(fftfreq[:-1],H_phase[refCH[0],pltidx,...])
    #     ax4.set_ylabel(r"$phase_{[%i,%i]}$"%(refCH[0],pltidx))
    #     ax4.set_xlabel("frequency in Hz")
    #     plt.show(block=False)

    # plt.show(block=True)
    for pltidx in range(nCHin):
        if pltidx == refCH[0]: continue
        plot1y = HdB[refCH[0],pltidx,:]
        plot2y = gamma2[refCH[0],pltidx,:]
        plot3y = IR.real[refCH[0],pltidx,:]
        plot4y = H_phase[refCH[0],pltidx,:]
        name1 = '$d, r \\text{ (solar radious)}$'
        freqs = fftfreq[:-1]
        tvec = tVec

    #     ax1.set_ylabel(r"$H1_{[%i,%i]}$ in dB"%(refCH[0],pltidx))
        trace1 = go.Scatter(
            x=freqs,
            y=plot1y,
            name='$H_{[%i,%i]}$ in dB'%(refCH[0],pltidx)
        )

        trace2 = go.Scatter(
            x=freqs,
            y=plot2y,
            name='$\\gamma^2_{[%i,%i]}$'%(refCH[0],pltidx)
        )

        trace3 = go.Scatter(
            x=tvec,
            y=plot3y,
            # name='$\\text{IR}_{[%i,%i]}$'%(refCH[0],pltidx)
            name=name1
        )

        trace4 = go.Scatter(
            x=freqs,
            y=plot4y,
            name='$Phase_{H,[%i,%i]}$'%(refCH[0],pltidx)
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
