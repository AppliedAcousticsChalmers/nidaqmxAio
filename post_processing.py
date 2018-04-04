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

def TFcalc(filename, blockSize, calibrationData):

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
    # Initial calculations
    keys = list(measurements.keys())
    dataCH = [i1 for i1, s1 in enumerate(keys) if "cDAQ" in s1]
    firstCH = dataCH[0]
    lastCH = dataCH[-1]
    keys = keys[firstCH:lastCH+1]
    refCH = [i2 for i2, s2 in enumerate(keys) if measurements['Reference_channel'] in s2][0]
    values = list(measurements.values())
    data = np.array(values[firstCH:lastCH+1])
    nCHin = int(len(data))
    print("Completed")

    # Removing the zero padding of the signal
    if simTime < 5 or "sweep" in signalType[0]:
        print("Removing the zero padding of the signal (This will take a while)... ", end="")
        convxymin = np.convolve((data[int(refCH), ...]), np.flipud(signalUnpadded), 'valid')
        convxymax = np.convolve(np.flipud(data[int(refCH), ...]), signalUnpadded, 'valid')

        pad_start = np.argmax(convxymin)
        pad_end = np.argmax(convxymax)
        data = data[..., pad_start:-pad_end]
        print("Completed")

    if "noise" in signalType:
        # Zeropadding
        numberOfsamples = int(len(data[0]))
        print("Zeropadding... ", end="")
        zeropad = blockSize - numberOfsamples % blockSize
        nBlocks = int(numberOfsamples // blockSize)
        data = np.pad(data, [[0, 0], [0, int(zeropad)]], 'constant', constant_values=0)
        dataMatrix = np.zeros((nCHin, blockSize, nBlocks))
        print("Completed")
        print("Signal padded with %i zeros at the end" % zeropad)


        # Arranging the data into blocks, Applying window and calculating the fft
        print("Arranging the data into blocks")
        for chidx in range(nCHin):
            i0 = 0
            pbar1 = tqdm(total=nBlocks)
            for idx in range(nBlocks):
                # Arranging into blocks
                dataMatrix[chidx, :, idx] = data[chidx, i0:i0+blockSize] #* win
                i0 += blockSize
                pbar1.update()
            pbar1.close()

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

        print("Saving the processed data... ", end="")
        fftfreq = np.fft.rfftfreq(blockSize, 1 / sr)
        tVec = np.linspace(0, len(IR[0, :]) / sr, len(IR[0, :]))
        measuredFRFs.update({"Hij": H, "HdBij": HdB, "HDij": HD, "IRij": IR, "gamma2ij": gamma2, "fftfreq": fftfreq, "tVec": tVec})
        print("Completed")

    elif "sweep" in signalType[0]:
        #Determining the start and stop frequencies of the signal
        f0 = int(signalType[1])
        f1 = int(signalType[2])
        print("Starting frequency f0 = %i Hz" % f0 )
        print("Ending frequency f1 = %i Hz" % f1 )

        #Main calculations
        print("Processing the data...", end="")
        H, HD, HdB, H_phase, IR, fftfreq, tVec = TF.deconvolution(data, sr, calibrationData, micAmp, refCH, f0, f1)
        print("Completed")

        #Saving the data
        print("Saving the processed data... ", end="")
        measuredFRFs.update({"Hij": H, "HdBij": HdB, "HDij": HD, "IRij": IR, "fftfreq": fftfreq, "tVec": tVec})
        print("Completed")

    for pltidx in range(nCHin):
        if pltidx == refCH:
            continue
        plot1y = HdB[pltidx, :]
        if "noise" in signalType: plot2y = gamma2[pltidx, :]
        plot3y = IR.real[pltidx, :]
        plot4y = H_phase[pltidx, :]
        freqs = fftfreq[:-1]
        tvec = tVec

    #     ax1.set_ylabel(r"$H1_{[%i,%i]}$ in dB"%(refCH[0],pltidx))
        trace1 = go.Scatter(
            x=freqs,
            y=plot1y,
            name='$H_{[%i,%i]}$ in dB' % (refCH, pltidx)
        )

        if "noise" in signalType:
           trace2 = go.Scatter(
               x=freqs,
               y=plot2y,
               name='$\\gamma^2_{[%i,%i]}$' % (refCH, pltidx)
           )

        trace3 = go.Scatter(
            x=tvec,
            y=plot3y,
            name='$\\text{IR}_{[%i,%i]}$' % (refCH, pltidx)
        )

        trace4 = go.Scatter(
            x=freqs,
            y=plot4y,
            name='$Phase_{H,[%i,%i]}$' % (refCH, pltidx)
        )

        fig = tools.make_subplots(rows=2, cols=2, subplot_titles=('HdB', 'gamma2', 'IR', 'H_phase'))

        fig.append_trace(trace1, 1, 1)
        if "noise" in signalType: fig.append_trace(trace2, 1, 2)
        fig.append_trace(trace3, 2, 1)
        fig.append_trace(trace4, 2, 2)

        fig['layout'].update(title=filename,
                             xaxis1=dict(title='frequency in Hz', type='log', autorange=True),
                             xaxis2=dict(title='frequency in Hz', type='log', autorange=True),
                             xaxis3=dict(title='time in s', type='lin', autorange=True),
                             xaxis4=dict(title='frequency in Hz', type='log', autorange=True)
                             )
        py.plot(fig, filename=filename + '.html')
    print("Process Finished")

    return measuredFRFs


def T60Shroeder(filename, blockSize, calibrationData):
    """Calculates the T_60 values in the bandwidth specified, using Schroeder's
       backwards intergration method.

       The function uses the calibration data, and the blockSize parameters to
       calculate the impulse responce of the room from the measured data folder,
       prior to the T_60 calculation.

       """
    # Calculating the IR
    print("Estimating the IR...", end="")
    measurements = TFcalc(filename, blockSize, calibrationData)
    print("Completed")

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

    print("Converting IR into 3rd octave bands")
    IR_3rd, fvec_3rd = fl.filters(IR,sr).band_filter(20,sr / 2.56,"third")
    IR_3rd = np.array(IR_3rd)[:,0,:]
    T_60 = []
    T_30 = []
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
        # plot1y = fl.amp2dB(IR_current)
        # plot2y = fl.amp2dB(envelope)
        # plot3y = fl.amp2dB(moav)
        # plot4y = fl.amp2dB(Rev_int)
        # plot5y = xr
        # plot6y = x_5dB
        # plot7y = x_35dB

        # trace1 = go.Scatter(
        #     x = tVec,
        #     y = plot1y,
        #     name='IR at %i Hz ' % (fvec_3rd[idx]),
        #     text = ['T_60 %i' % (T_60[idx])],
        #     textposition='bottom'
        # )

        # trace2 = go.Scatter(
        #     x=tVec,
        #     y=plot2y,
        #     name='envelope %i Hz' % (fvec_3rd[idx])
        # )

        # trace3 = go.Scatter(
        #     x=tVec,
        #     y=plot3y,
        #     name='moav %i Hz' % (fvec_3rd[idx])
        # )

        # trace4 = go.Scatter(
        #     x=tVec,
        #     y=plot4y,
        #     name='Rev intH %i Hz' % (fvec_3rd[idx])
        # )

        # trace5 = go.Scatter(
        #     x=tVec,
        #     y=plot5y,
        #     name='fit %i Hz' % (fvec_3rd[idx])
        # )

        # trace6 = go.Scatter(
        #     x=tVec,
        #     y=plot6y,
        #     name='-5dB'
        # )

        # trace7 = go.Scatter(
        #     x=tVec,
        #     y=plot7y,
        #     name='-35dB'
        # )
        # data = [trace1, trace2, trace4]#, trace3, trace5, trace6, trace7]
        # py.plot(data, filename=filename + '.html')

    #Saving the data
    RT = {'T60Schroeder': T_60, 'T30Schroeder': T_30, 'fVec_3rd': fvec_3rd}

    return [measurements, RT]
