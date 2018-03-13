#Python libraries
# import matplotlib
# import matplotlib.pyplot as plt
from pyqtgraph.Qt import QtGui
import pyqtgraph as pg
import numpy as np
import time
from tqdm import tqdm
#Program libraries
import utils as gz
import meas_signal as sig
import post_processing as pp
#nidaqmx libraries
import collections
import nidaqmx.system
from nidaqmx.task import Task
from nidaqmx.constants import AcquisitionType, TaskMode, Coupling
from nidaqmx.stream_writers import AnalogMultiChannelWriter
from nidaqmx.stream_readers import AnalogMultiChannelReader
system = nidaqmx.system.System.local()


#Data acquisition script
def ni_io_tf(args):
    #Getting the initial arguments
    sim_time = args.time
    sample_rate = args.sampleRate
    if sample_rate % 2 !=0:
        sample_rate -= sample_rate % 2
        sample_rate = int(sample_rate)
    bufferSize = int(sample_rate // 2)
    signal_temp = sig.create_signal(args.sampleRate, args.time, args.pad_samples, args.signalType, args.aoRange)
    signal = signal_temp[0]
    signal_unpadded = signal_temp[1]
    number_of_channels_in = args.channelsIn
    number_of_channels_out = args.channelsOut
    ai_range = args.aiRange
    ao_range = args.aoRange
    calibrationData = args.calData
    number_of_samples = sim_time*sample_rate + args.pad_samples
    number_of_samples += bufferSize - (number_of_samples%bufferSize)
    if not np.sum(number_of_channels_out)<= 1: signal = np.tile(signal, [np.sum(number_of_channels_out), 1])
    micAmp = input("Please specify the mic preamplification (eg: 1, 10, 100, default = 1)")
    if not micAmp: micAmp = 1
    else: micAmp = int(micAmp)
    #Creating the dictionary to store the data
    measurements = {'simulationTime': sim_time, 'SampleRate': sample_rate, 'Signal': args.signalType, 'Input_range_in_Vrms': ai_range, 'Output_range_in_Vrms': ao_range, 'micAmp': micAmp, 'Unpadded_signal': signal_unpadded, 'Reference_channel':[]}
    measuredFRFs = measurements

    #Reading the pressent in/out channels
    channel_list = []
    for device in system.devices:
        channel_list.append({"ai":[channels_ai.name for _, channels_ai in enumerate(device.ai_physical_chans)],
                            "ao":[channels_ao.name for _, channels_ao in enumerate(device.ao_physical_chans)]})

    #Selecting the channels that were selected by the user
    idx_ai = []
    idx_ao = []
    if np.sum(number_of_channels_in) == 1:
        values_read = np.empty((0))
    else:
        values_read = np.empty((np.sum(number_of_channels_in),0))
    for idx in range(len(channel_list)):
        if channel_list[idx]['ai']!=[] : idx_ai.append(idx)
        if channel_list[idx]['ao']!=[] : idx_ao.append(idx)
    idx_ai = idx_ai[:len(number_of_channels_in)]
    if sum(number_of_channels_in)==0: idx_ai=[]
    idx_ao = idx_ao[:len(number_of_channels_out)]
    if sum(number_of_channels_out)==0: idx_ao=[]

    # Selecting the input channel that is going to be used as reference
    ch_count = 0
    chSelect = []
    print("Please select which channel should be used as reference in post processing. (pressing enter will default to the first channel in the list)")
    for idx, device_idx in enumerate(idx_ai):
        if number_of_channels_in[idx] == 0 or idx_ai==[] : continue
        for ch_num in range(number_of_channels_in[idx]):
            chSelect.append(channel_list[device_idx]['ai'][ch_num])
            print("[" + str(ch_count)+ "]"+" "+chSelect[ch_count])
            ch_count += 1
    selection = input()
    if selection: selection = int(selection); refChannel = chSelect[selection]
    else: refChannel = chSelect[0]; selection = 0
    measurements.update({'Reference_channel':refChannel})

    #Setting up the in/out task
    Coupling.AC
    with nidaqmx.Task() as write_task, nidaqmx.Task() as read_task:
        for idx, device_idx in enumerate(idx_ao):
            if number_of_channels_out[idx] == 0 or idx_ao==[] : continue
            write_task.ao_channels.add_ao_voltage_chan(channel_list[device_idx]['ao'][0]+":%i"%(number_of_channels_out[idx]-1), max_val=ao_range, min_val=-ao_range)

        for idx, device_idx in enumerate(idx_ai):
            if number_of_channels_in[idx] == 0 or idx_ai==[] : continue
            read_task.ai_channels.add_ai_voltage_chan(channel_list[device_idx]['ai'][0]+":%i"%(number_of_channels_in[idx]-1), max_val=ai_range, min_val=-ai_range)
        if idx_ai != []:
            read_task.timing.cfg_samp_clk_timing(rate=sample_rate, sample_mode = AcquisitionType.CONTINUOUS)

        if idx_ao != []:
            write_task.timing.cfg_samp_clk_timing(rate=sample_rate, sample_mode = AcquisitionType.CONTINUOUS)
            write_task.write(signal)
        if idx_ao: write_task.control(TaskMode.TASK_COMMIT)
        if idx_ai: read_task.control(TaskMode.TASK_COMMIT)
        if idx_ao: write_task.start()
        if idx_ai: read_task.start()

        #Starting the data aquition/reproduction

        #Intialization
        tVec = np.linspace(0, int(bufferSize / sample_rate), int(bufferSize))
        fftfreq = np.fft.rfftfreq(bufferSize, 1 / sample_rate)
        timeCounter = 0
        blockidx = 0
        # previous_buffer = []
        # current_buffer = []
        values_read = np.zeros((sum(number_of_channels_in),int(number_of_samples)))
        current_buffer = np.zeros((sum(number_of_channels_in),int(bufferSize)))
        previous_buffer = np.zeros((sum(number_of_channels_in),int(bufferSize)))

        #Figure creation
        from pyqtgraph.Qt import QtGui
        app = QtGui.QApplication([])
        import pyqtgraph as pg

        win = pg.GraphicsLayoutWidget()
        win.setWindowTitle('And awaaaaay we go!')
        win.resize(1000, 600)
        win.show()

        pg.setConfigOptions(antialias=False)
        p = []
        curve = []
        plotCounter = 0
        if sum(number_of_channels_in) == 1:
            numPlots = 1
            plotsPerRow = 1
        elif sum(number_of_channels_in) == 2:
            numPlots = 4
            plotsPerRow = 2
        else:
            numPlots =  sum(number_of_channels_in ) * 4
            plotsPerRow = 4
        downsample = numPlots*0+1
        for i in range(numPlots):
            # if i==selection and sum(number_of_channels_in) !=1:j = 1; continue
            # else: j = i
            if plotCounter == 0:
                p.append(win.addPlot(title='Spectum Ch:'+str(i)))
                p[i].showGrid(True, True)
                # p[i].setRange(yRange=[-100, 0])
                if sum(number_of_channels_in) > 1:
                    p[i].setLogMode(True, False)
                    curve.append(p[i].plot(fftfreq, np.zeros(len(fftfreq)), pen=(173,255,47)))
                    curve.append(p[i].plot(fftfreq, np.zeros(len(fftfreq)), pen=(200,200,200)))
                else:
                    curve.append(p[i].plot(tVec, np.zeros(len(tVec)), pen=(200,200,200)))
                plotCounter += 1
            elif plotCounter == 1:
                p.append(win.addPlot(title='H Ch: '+str(i)))
                p[i].setLogMode(True, False)
                p[i].showGrid(True, True)
                # p[i].setRange(yRange=[-100, 0])
                curve.append(p[i].plot(fftfreq, np.zeros(len(fftfreq)), pen=(200, 200, 200)))
                plotCounter += 1
                if sum(number_of_channels_in) == 2: win.nextRow();
            elif plotCounter == 2:
                p.append(win.addPlot(title='IR Ch: '+str(i)))
                p[i].setLogMode(False, False)
                p[i].showGrid(True, True)
                # p[i].setRange(yRange=[-0.05, 0.05])
                curve.append(p[i].plot(tVec, np.zeros(len(tVec)), pen=(200, 200, 200)))
                plotCounter += 1
            else:
                p.append(win.addPlot(title='gamma2 Ch: '+str(i)))
                p[i].setLogMode(True, False)
                p[i].showGrid(True, True)
                p[i].setRange(yRange=[0, 1.1])
                curve.append(p[i].plot(fftfreq, np.zeros(len(fftfreq)), pen=(200, 200, 200)))
                plotCounter += 1
            if (i+1) % plotsPerRow == 0 and sum(number_of_channels_in) !=2:
                win.nextRow()
                plotCounter = 0

        #Main loop
        if idx_ai:
            pbar_ai = tqdm(total=number_of_samples)
            # if number_of_channels_in[0] >= 2:
            while timeCounter < number_of_samples:
                current_buffer = read_task.read(number_of_samples_per_channel=bufferSize)
                current_buffer = np.array(current_buffer)
                values_read[:,timeCounter:timeCounter+bufferSize] =  current_buffer
                if number_of_channels_in[0] >= 2:
                    previous_buffer, spectra, blockidx, H, HdB, H_phase, IR, gamma2  = pp.h1_estimator_live(current_buffer, previous_buffer, blockidx, calibrationData, micAmp, selection)
                else:spectra = np.array(current_buffer)
                    #Plotting
                for i in range(0,numPlots,4):
                    if numPlots > 1:
                        curve[i].setData(gz.amp2db(spectra[0,0:int(bufferSize//2)]), antialias=True, downsample=downsample, downsampleMethod='subsample')
                        curve[i+1].setData(gz.amp2db(spectra[1,0:int(bufferSize//2)]), antialias=True, downsample=downsample, downsampleMethod='subsample')
                        curve[i+2].setData(HdB[1,...], antialias=True, downsample=downsample, downsampleMethod='subsample')
                        curve[i+3].setData(IR.real[1,...], antialias=True, downsample=downsample, downsampleMethod='subsample')
                        curve[i+4].setData(gamma2[1,...], antialias=True, downsample=downsample, downsampleMethod='subsample')
                    else:
                        curve[i].setData(spectra, antialias=True, downsample=downsample, downsampleMethod='subsample')
                    pg.QtGui.QApplication.processEvents()

                timeCounter += bufferSize
                pbar_ai.update(bufferSize)
            pbar_ai.close()
        elif idx_ao:
            time.sleep(sim_time)
            pbar_ao = tqdm(total=sim_time)
            for pbar_ao_idx in range(pbar_ao):
                pbar_ao.update()
            pbar_ao.close()

    #Updating the dictionary with the measured data
    ch_count = 0
    for idx, device_idx in enumerate(idx_ai):
        if number_of_channels_in[idx] == 0 or idx_ai==[] : continue
        for ch_num in range(number_of_channels_in[idx]):
            measurements.update({(channel_list[device_idx]['ai'][ch_num]):np.array(values_read)[ch_count, ...]})
            ch_count += 1

    if number_of_channels_in[0] >= 2: measuredFRFs.update({"Hij":H, "HdBij":HdB, "H_phase_ij": H_phase, "IRij":IR, "gamma2ij":gamma2, "fftfreq":fftfreq, "tVec":tVec})


    return measurements
