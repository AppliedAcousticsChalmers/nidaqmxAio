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
import roomAcousticsModule as RA
import filtersAndMathtools as flm
import plotting as myplt


def mic_calibration(measurements, sensitivity):
    # Reading the configuration
    keys = list(measurements.keys())
    sr = measurements['SampleRate']
    sensitivity = float(sensitivity) * float(1e-3)
    calMeasurement = min([i1 for i1, s1 in enumerate(keys) if "cDAQ" in s1])
    values = list(measurements.values())
    data = np.array(values[calMeasurement]) / sensitivity

    #
    Lpmeas = 20 * np.log10( np.sqrt( np.mean( data ** 2 )) / ( 2e-5  ))
    calibration = [10 ** ((93.9791 - Lpmeas)/20), sensitivity, {"data":data,"SampleRate": sr}]
    #

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

    # # Removing the zero padding of the signal
    # if simTime < 5 or ("sweep" in signalType[0]):
    #     print("Removing the zero padding of the signal (This will take a while)... ", end="")
    #     convxymin = np.convolve((data[int(refCH), ...]), np.flipud(signalUnpadded), 'valid')
    #     convxymax = np.convolve(np.flipud(data[int(refCH), ...]), signalUnpadded, 'valid')

    #     pad_start = np.argmax(convxymin)
    #     pad_end = np.argmax(convxymax)
    #     data = data[..., pad_start:]
    #     data = data[..., pad_start:-pad_end]
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
            previous_buffer, _, blockidx, H, HD, HdB, H_phase, IR, gamma2, SNR = TF.h1_estimator(current_buffer, previous_buffer, blockidx, calibrationData, micAmp, refCH)
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
                plots = {'TF':[fftfreq, HdB[pltidx, :]],
                         'gamma2':[fftfreq, gamma2[pltidx, :]],
                         'IR':[tVec, IR[pltidx, :]],
                         'phase':[fftfreq, H_phase[pltidx, :]]}
            elif ("sweep" in signalType[0]):
                plots = {'TF':[fftfreq, HdB[pltidx, :]],
                         'spectrogram':[spectrogramm[pltidx][1], spectrogramm[pltidx][0],spectrogramm[pltidx][2]],
                         'IR':[tVec, IR.real[pltidx, :]],
                         'phase':[fftfreq, H_phase[pltidx, :]]}

            myplt.irSubPlot(plots, filename + "_" + str(pltidx), "Channel: " + keys[pltidx])

    print("Process Finished")
    print("")

    return measuredFRFs


def RAC(filename, blockSize, f_min, f_max, bandWidth, calibrationData, plotting):
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


    db = lambda x: flm.NormAmp2dB(x)

    T_30, T_60, tau, fVec_band, IRs, IR2s, Rev_int_curves, line_fits_T30, line_fits_T60 = RA.T60(IR[0], sr, tVec, f_min=float(f_min), f_max=float(f_max), bandWidth=str(bandWidth))

    # #Plotting

    for idx in range(0,len(T_60)):
        # Scroeder one band
        plot1y = db(IRs[idx])
        plot2y = db(IR2s[idx])
        plot3y = db(Rev_int_curves[idx])
        plot4y = line_fits_T30[idx][0]
        plot5y = line_fits_T60[idx][0]
        plot6y = line_fits_T30[idx][1]
        plot7y = line_fits_T30[idx][2]
        plot8y = line_fits_T60[idx][2]
        plots2 = {'0':[tVec,plot1y,'Impulse Response at %i Hz 3rd octave band: ' % (fVec_band[idx]) ],
            '1':[tVec,plot2y,'Squared IR'],
            '2':[tVec,plot3y,'Schroeder curve'],
            '3':[tVec,plot4y,'Linear fit T30'],
            '4':[tVec,plot5y,'Linear fit T60'],
            '5':[tVec,plot6y,'-5dB'],
            '6':[tVec,plot7y,'-35dB'],
            '7':[tVec,plot8y,'-65dB'],
             'xAxisTitle': 'time in s',
             'yAxisTitle': 'Amplitude',
             'plotTitle': filename + "_@_" + str(fVec_band[idx]),
             'scale':'lin'}

        if 'T60_one_band' in plotting: myplt.singlePlot(plots2)

    # T60 and T30 in bands
    plots2 = {'0': [fVec_band,T_60,"T60"],
              '1': [fVec_band,T_30,"T30"],
              'xAxisTitle': 'frequency in Hz',
              'yAxisTitle': 'time in s',
              'plotTitle': filename + "_RT",
              'scale':'log'}

    if 'T60_3rd' in plotting: myplt.singlePlot(plots2)

    # Saving the data
    print("Saving the processed data... ", end="")
    RT = {'T60Schroeder': T_60, 'T30Schroeder': T_30, 'fVec_band': fVec_band,
          "tVec":tVec, 'fVec':fVec, 'IRs':IRs, "IR2s":IR2s, "Rev_int_curves":Rev_int_curves, "line_fits_T30":line_fits_T30, "line_fits_T60":line_fits_T60}
    print("Completed")

    print("Process Finished")
    print("")

    return [measurements, RT]
