# Python libraries
from plotly import tools
import plotly.offline as py
import plotly.graph_objs as go
import numpy as np
import math
from scipy import signal, array, polyfit, polyval
from scipy.signal import hilbert
from tqdm import tqdm
# Program libraries
import utils as gz
import TFestimation as TF
import filters as fl
import plotting as myplt


def mic_calibration(measurements, sensitivity):
    # Reading the configuration
    keys = list(measurements.keys())
    sensitivity = float(sensitivity) * float(1e-3)
    calMeasurement = min([i1 for i1, s1 in enumerate(keys) if "cDAQ" in s1])
    values = list(measurements.values())
    data = np.array(values[calMeasurement])
    Lpmeas = max(gz.amp2db((data)/(sensitivity), ref=2e-5))
    calibration = [10 ** ((93.9794 - Lpmeas)/20), sensitivity]

    return calibration

def TFcalc(filename, blockSize, calibrationData, plotting):
    print("Transfer function estimation")

    # Reading the configuration
    print("Importing raw data... ", end="")
    inputfile = np.load(filename)
    measurements = inputfile.item()
    measuredFRFs = measurements
    simTime = measurements["simulationTime"]
    sr = measurements["SampleRate"]
    micAmp = measurements["micAmp"]
    blockSize = int(blockSize)  # is gonna be handled as a command line input
    signalUnpadded = measurements["Unpadded_signal"]
    signalType = measurements["Signal"]

    # Locate the part of the file containing the measured data
    keys = list(measurements.keys())
    dataCH = [i1 for i1, s1 in enumerate(keys) if "cDAQ" in s1]
    firstCH = dataCH[0]
    lastCH = dataCH[-1]
    keys = keys[firstCH:lastCH+1]
    # Locate the reference channel or ask for user input if none is given
    if measurements['Reference_channel'] == 'none':
        print("Please select which channel should be used as reference in post processing. (pressing enter will default to the first channel in the list)")
        idx = 0
        for name in keys:
            print('[' + str(idx) + '] ' + name)
            idx += 1
        selection = input()
        if not selection: selection = 0
        measurements['Reference_channel'] = keys[int(selection)]
    refCH = [i2 for i2, s2 in enumerate(keys) if measurements['Reference_channel'] in s2][0]

    # Store the measured data into an array and count the total number of channels
    values = list(measurements.values())
    data = np.array(values[firstCH:lastCH+1])
    nCHin = int(len(data))
    print("Completed")

    # Removing the zero padding of the signal
    # if simTime < 5 or ("sweep" in signalType[0]):
        # print("Removing the zero padding of the signal (This will take a while)... ", end="")
        # convxymin = np.convolve((data[int(refCH), ...]), np.flipud(signalUnpadded), 'valid')
        # convxymax = np.convolve(np.flipud(data[int(refCH), ...]), signalUnpadded, 'valid')

        # pad_start = np.argmax(convxymin)
        # pad_end = np.argmax(convxymax)
        # data = data[..., pad_start:]
        # data = data[..., pad_start:-pad_end]
        # print("Completed")


    # Calculating the TF for Noise signals
    if ("noise" in signalType[0]):

        # Zeropadding
        numberOfsamples = int(len(data[0]))
        print("Zeropadding... ", end="")
        zeropad = blockSize - numberOfsamples % blockSize
        nBlocks = int(numberOfsamples // blockSize)
        data = np.pad(data, [[0, 0], [0, int(zeropad)]], 'constant', constant_values=0)
        dataMatrix = np.zeros((nCHin, blockSize, nBlocks))
        print("Completed")
        print("Signal padded with %i zeros at the end" % zeropad)

        # Arranging the data into blocks
        print("Arranging the data into blocks")
        for chidx in range(nCHin):
            print(np.shape(data))
            i0 = 0
            pbar1 = tqdm(total=nBlocks)
            for idx in range(nBlocks):
                # Arranging into blocks
                dataMatrix[chidx, :, idx] = data[chidx, i0:i0+blockSize]
                i0 += blockSize
                pbar1.update()
            pbar1.close()

        # Calculations using the H1 estimator
        blockidx = 0
        current_buffer = np.zeros((nCHin, int(blockSize)))
        previous_buffer = np.zeros((nCHin, int(blockSize)))
        print("Processing the data")
        pbar2 = tqdm(total=nBlocks)
        for idx in range(nBlocks):
            current_buffer = dataMatrix[:,:,idx]
            previous_buffer, _, blockidx, H, HD, HdB, H_phase, IR, gamma2 = TF.h1_estimator(current_buffer, previous_buffer, blockidx, calibrationData, micAmp, refCH)
            pbar2.update()
        pbar2.close()

        # Saving the data
        print("Saving the processed data... ", end="")
        fftfreq = np.fft.rfftfreq(blockSize, 1 / sr)
        tVec = np.linspace(0, len(IR[0, :]) / sr, len(IR[0, :]))
        measuredFRFs.update({"Hij": H, "HdBij": HdB, "HDij": HD, "IRij": IR, "gamma2ij": gamma2, "fftfreq": fftfreq, "tVec": tVec})
        print("Completed")

    # Calculating the TF for Sweep signal
    elif ("sweep" in signalType[0]):

        # Determining the start and stop frequencies of the signal
        f0 = int(signalType[1])
        f1 = int(signalType[2])
        print("Starting frequency f0 = %i Hz" % f0 )
        print("Ending frequency f1 = %i Hz" % f1 )

        #Main calculations
        print("Processing the data...", end="")
        H, HD, HdB, H_phase, IR, fftfreq, tVec, spectrogramm = TF.deconvolution(data, sr, blockSize, calibrationData, micAmp, refCH, f0, f1)
        print("Completed")

        # Saving the data
        print("Saving the processed data... ", end="")
        measuredFRFs.update({"Hij": H, "HdBij": HdB, "HDij": HD, "IRij": IR, "H_phase":H_phase, "fftfreq": fftfreq, "tVec": tVec})
        print("Completed")


    # Plotting
    if 'TF' in plotting:
        for pltidx in range(nCHin):
            if pltidx == refCH:
                continue
            if ("noise" in signalType[0]):
                plots = {'TF':[fftfreq[:-1], HdB[pltidx, :]],
                         'gamma2':[fftfreq[-1], gamma2[pltidx, :]],
                         'IR':[tVec, IR.real[pltidx, :]],
                         'phase':[fftfreq[-1], H_phase[pltidx, :]]}
            elif ("sweep" in signalType[0]):
                plots = {'TF':[fftfreq[:-1], HdB[pltidx, :]],
                         'spectrogram':[spectrogramm[pltidx][1], spectrogramm[pltidx][0],spectrogramm[pltidx][2]],
                         'IR':[tVec, IR.real[pltidx, :]],
                         'phase':[fftfreq[-1], H_phase[pltidx, :]]}

            myplt.irSubPlot(plots, filename + "_" + str(pltidx), "Channel: " + keys[pltidx])

    # Time signals plotting
    if 'timeSig' in plotting:
        time_signals_plot = {'xAxisTitle': 'time in s',
                             'yAxisTitle': 'Amplitude',
                             'plotTitle':  filename + '_time_signals',
                             'scale':'lin'}

        full_tVec = np.linspace(0,simTime, int(len(data[0])))
        time_signals = data
        time_signals[:refCH, :] = (data[:refCH, :] / micAmp) * calibrationData[0]
        time_signals[refCH+1:, :] = (data[refCH+1:, :] / micAmp) * calibrationData[0]

        for idx in range(nCHin):
            time_signals_plot.update({str(idx):[full_tVec, time_signals[idx,:],"Channel:" + keys[idx]]})

        myplt.singlePlot(time_signals_plot)

    print("Process Finished")

    return measuredFRFs


def T60Shroeder(filename, blockSize, calibrationData, plotting):
    """Calculates the T_60 values in the bandwidth specified, using Schroeder's
       backwards intergration method.

       The function uses the calibration data, and the blockSize parameters to
       calculate the impulse responce of the room from the measured data folder,
       prior to the T_60 calculation.

       """
    # Calculating the IR
    measurements = TFcalc(filename, blockSize, calibrationData, plotting)

    print("T60 Estimation")
    #Importing the required data
    print("Importing the required data... ", end="")
    sr = measurements["SampleRate"]
    IR = measurements["IRij"]
    tVec = measurements["tVec"]
    fVec = measurements["fftfreq"]
    print("Completed")

    # Reading the configuration
    print("Reading the configuration...", end="")
    keys = list(measurements.keys())
    dataCH = [i1 for i1, s1 in enumerate(keys) if "cDAQ" in s1]
    firstCH = dataCH[0]
    lastCH = dataCH[-1]
    keys = keys[firstCH:lastCH+1]
    refCH = [i2 for i2, s2 in enumerate(keys) if measurements['Reference_channel'] in s2][0]
    IR = np.delete(IR, (refCH), axis=0)
    print("completed")

    #Removing the source - rec distance
    print("Removing source - receiver distance...", end="")
    # for idx in range(0, len(IR)):
    IR_head, IR_tail = fl.filters(IR[0,:], sr).removeDelay()

    IR_headVEC = np.array(range(0,len(IR_head)))
    IR_tailVEC = len(IR_head) + np.array(range(0,len(IR_tail)))
    print("completed")


    print("Converting IR into 3rd octave bands...", end="")
    IR_3rd, fvec_3rd = fl.filters(IR_tail,sr).band_filter(20,sr / 2.56,"third")
    IR_3rd = np.array(IR_3rd)
    print("completed")

    T_60 = []
    T_30 = []
    print("Calculating T60 for each band")
    pbar1 = tqdm(total=len(IR_3rd-1))
    for idx in range(0,len(IR_3rd-1)):

        IR_current = IR_3rd.real[idx,:]

        #Finding the analytical signal and taking it's envelope through the hilbert transform
        s_anal = hilbert(IR_current)
        envelope =  abs(s_anal)

        #Calculating the moving average of the envelope
        moav = fl.smooth(envelope, int(0.1 * sr), 'flat')
        moav_dB = fl.amp2dB(moav)

        #Finding the correct threshold for the calculation of tau
        threshold = np.average(moav_dB[int(len(moav)*0.75):int(len(moav_dB)*0.75)+2000])
        tau = fl.findTau(moav_dB, tVec, threshold)

        #Calculating the Shroeder curve
        Rev_int = np.zeros((len(IR_current[:tau])))
        Rev_int[tau::-1] =(np.cumsum(moav[tau:0:-1]) / np.sum(moav[0:tau]))
        Rev_int_dB = fl.amp2dB(Rev_int)

        #Fitting a line in the Shroeder curve
        fit_range = int(len(Rev_int_dB) * 0.75)
        x = tVec[:fit_range]
        Rev_int_fit = Rev_int_dB[:fit_range]
        (b, a) = polyfit(x, Rev_int_fit, 1)
        xr=polyval([b,a],tVec)

        #Fitting the guiding lines
        x_5dB = polyval([0,max(Rev_int_dB)-5],tVec)
        x_35dB = polyval([0,max(Rev_int_dB)-35],tVec)

        #Calculating T_60
        T_60.append( -60 / b)

        #Calculating T_30
        t_5dB = tVec[np.argmin(abs(xr-x_5dB))]
        t_35dB = tVec[np.argmin(abs(xr-x_35dB))]
        T_30.append(t_35dB - t_5dB)



    # #Plotting

        # Scroeder one band
        plot1y = fl.amp2dB(IR_current)
        plot2y = fl.amp2dB(envelope)
        plot3y = fl.amp2dB(moav)
        plot4y = fl.amp2dB(Rev_int)
        plot5y = xr
        plot6y = x_5dB
        plot7y = x_35dB
        plots2 = {'0':[tVec,plot1y,'Impulse Responce at %i Hz 3rd octave band: ' % (fvec_3rd[idx]) ],
             '1':[tVec,plot2y,'Envelope'],
             '2':[tVec,plot3y,'Moving average filter'],
             '3':[tVec,plot4y,'Schroeder curver'],
             '4':[tVec,plot5y,'Linear fit'],
             '5':[tVec,plot6y,'-5dB'],
             '6':[tVec,plot7y,'-35dB'],
             'xAxisTitle': 'time in s',
             'yAxisTitle': 'Amplitude',
             'plotTitle': filename + "_@_" + str(fvec_3rd[idx]),
             'scale':'lin'}

        if 'T60_one_band' in plotting: myplt.singlePlot(plots2)

        pbar1.update()
    pbar1.close()

    # Source - Receiver removal
    plots = {'0':[np.array(len(IR[0,:])),IR.real[0,:],'Impulse Responce'],
             '1':[IR_headVEC,IR_head.real,'head'],
             '2':[IR_tailVEC,IR_tail.real,'tail'],
             'xAxisTitle': 'time in s',
             'yAxisTitle': 'Amplitude',
             'plotTitle': 'initial delay removal',
             'scale':'lin'}

    if 'T60_SRR' in plotting: myplt.singlePlot(plots)

    # T60 final
    plots3 = {'0': [fvec_3rd,T_60,"T60"],
              'xAxisTitle': 'frequency in Hz',
              'yAxisTitle': 'time in s',
              'plotTitle': filename + "T_60",
              'scale':'log'}

    if 'T60_3rd' in plotting: myplt.singlePlot(plots3)
    #Saving the data
    RT = {'T60Schroeder': T_60, 'T30Schroeder': T_30, 'fVec_3rd': fvec_3rd}

    return [measurements, RT]
