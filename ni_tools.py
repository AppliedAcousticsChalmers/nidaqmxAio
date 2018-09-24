# Python libraries
from pyqtgraph.Qt import QtGui
import pyqtgraph as pg
import numpy as np
import time
from tqdm import tqdm
import copy
# Program libraries
import systemUtils as syU
import utils as gz
import meas_signal as sig
import TFestimation as TF
import plotting as myPlt
# nidaqmx libraries
import nidaqmx.system
import nidaqmx
from nidaqmx.task import Task
from nidaqmx.constants import AcquisitionType, TaskMode, Coupling
from nidaqmx.stream_writers import AnalogMultiChannelWriter
from nidaqmx.stream_readers import AnalogMultiChannelReader
system = nidaqmx.system.System.local()

# Data acquisition script
def ni_io_tf(args, calibrationData=[1, 1], cal=False):
    if cal == True:
        input_sim_time = 10
        input_sample_rate = 51200
        input_cutoffTime = 0
        bufferSize = 5120
        number_of_channels_in = [1]
        number_of_channels_out = [0]
        pad_samples = 0
    else:
        bufferSize = args.bufferSize
        input_sim_time = args.time
        input_sample_rate = args.sampleRate
        input_cutoffTime = args.cutoffTime
        bufferSize = args.bufferSize
        number_of_channels_in = args.channelsIn
        number_of_channels_out = args.channelsOut
        pad_samples = args.pad_samples

    # Setting the sampling frequency to be an even number, so the Niquist frequency is at the last bin of each block
    # Correcting the sampling rate to correspond to one supported by the cDAQ system's clock (see sampleRateCorrection funtion for more information)
    if input_sample_rate % 2 != 0 and cal == False:
        sample_rate -= input_sample_rate % 2
        sample_rate = int(syU.sampleRateCorrection(sample_rate))
    else:
        sample_rate = int(syU.sampleRateCorrection(input_sample_rate))
    if input_sample_rate != sample_rate:
        print("Given sample rate of:",input_sample_rate,"Hz is not supported by the Hardware.")
        print("Switching to:",sample_rate,"Hz")

    #Signal creation
    signal_type = args.signalType[0].split("_")
    signal_temp = sig.create_signal(sample_rate, input_sim_time, pad_samples, args.signalType, args.aoRange, input_cutoffTime)
    signal = signal_temp[0]
    signal_unpadded = signal_temp[1]
    if not np.sum(number_of_channels_out) <= 1:
        signal = np.tile(signal, [np.sum(number_of_channels_out), 1])
    number_of_samples = signal_temp[2]
    lastBufferZeroPadding = bufferSize - (number_of_samples % bufferSize)
    sim_time = signal_temp[3]
    cutoffTime = signal_temp[4]
    if sim_time != input_sim_time or cutoffTime != input_cutoffTime:
        print("In order to avoid rounding errors the following values were adjusted as:")
        if sim_time != input_sim_time:print("Simulation time used:", sim_time)
        if cutoffTime != input_cutoffTime: print("cutoffTime time used:", cutoffTime)

    # Reading the in/out range
    ai_range = args.aiRange
    ao_range = args.aoRange

    # Checking the microphone amplification
    micAmp = args.micAmp
    if not micAmp: micAmp = 1
    else: micAmp = int(micAmp)

    # Creating the dictionary to store the data
    measurements = {'simulationTime': sim_time, 'SampleRate': float(sample_rate), 'Signal': args.signalType, 'Input_range_in_Vrms': ai_range, 'Output_range_in_Vrms': ao_range, 'bufferSize': bufferSize, 'micAmp': micAmp, 'Unpadded_signal': signal_unpadded, 'Reference_channel': [], 'Comments': args.comment}

    # Reading the present in/out channels
    channel_list = []
    for device in system.devices:
        channel_list.append({"ai": [channels_ai.name for _, channels_ai in enumerate(device.ai_physical_chans)],
                             "ao":  [channels_ao.name for _, channels_ao in enumerate(device.ao_physical_chans)]})

    # Channels selected by user
    idx_ai = []
    idx_ao = []
    if np.sum(number_of_channels_in) == 1:
        values_read = np.empty((0))
    else:
        values_read = np.empty((np.sum(number_of_channels_in), 0))
    for idx in range(len(channel_list)):
        if channel_list[idx]['ai'] != []:
            idx_ai.append(idx)
        if channel_list[idx]['ao'] != []:
            idx_ao.append(idx)
    idx_ai = idx_ai[:len(number_of_channels_in)]
    if sum(number_of_channels_in) == 0:
        idx_ai = []
    idx_ao = idx_ao[:len(number_of_channels_out)]
    if sum(number_of_channels_out) == 0:
        idx_ao = []

    # Reference input channel selection
    ch_count = 0
    chSelect = []
    if idx_ai != []:
        if sum(number_of_channels_in) > 1 and not args.refferenceChannel:print("Please select which channel should be used as reference in post processing. (pressing enter will default to the first channel in the list)")
        for idx, device_idx in enumerate(idx_ai):
            if number_of_channels_in[idx] == 0 or idx_ai == []:
                continue
            for ch_num in range(number_of_channels_in[idx]):
                chSelect.append(channel_list[device_idx]['ai'][ch_num])
                print("[" + str(ch_count) + "]" + " " + chSelect[ch_count])
                ch_count += 1
        if sum(number_of_channels_in) > 1:
            print("[" + str(ch_count) + "]" + " " + "No reference channel")
            if not args.refferenceChannel:selection = input()
            else: selection = args.refferenceChannel
            if selection:
                if int(selection) == int(ch_count):
                    selection = "none"
                    refChannel = "none"
                else:
                    selection = int(selection)
                    refChannel = chSelect[selection]
            else:
                refChannel = chSelect[0]
                selection = 0
        else:
            selection = "none"
            refChannel = "none"
        measurements.update({'Reference_channel': refChannel})

    # Setting up the in/out task
    Coupling.AC
    with nidaqmx.Task() as write_task, nidaqmx.Task() as read_task:
        for idx, device_idx in enumerate(idx_ao):
            if number_of_channels_out[idx] == 0 or idx_ao == []:
                continue
            write_task.ao_channels.add_ao_voltage_chan(channel_list[device_idx]['ao'][0]+":%i" % (number_of_channels_out[idx]-1), max_val=ao_range, min_val=-ao_range)

        for idx, device_idx in enumerate(idx_ai):
            if number_of_channels_in[idx] == 0 or idx_ai == []:
                continue
            read_task.ai_channels.add_ai_voltage_chan(channel_list[device_idx]['ai'][0]+":%i" % (number_of_channels_in[idx]-1), max_val=ai_range, min_val=-ai_range)
        if idx_ai != []:
            read_task.timing.cfg_samp_clk_timing(rate=sample_rate, sample_mode=AcquisitionType.CONTINUOUS)

        if idx_ao != []:
            write_task.out_stream.regen_mode = nidaqmx.constants.RegenerationMode.DONT_ALLOW_REGENERATION
            write_task.timing.cfg_samp_clk_timing(rate=sample_rate, sample_mode=AcquisitionType.CONTINUOUS)
            write_task.write(signal)

        if idx_ao:
            write_task.control(TaskMode.TASK_COMMIT)
        if idx_ai:
            read_task.control(TaskMode.TASK_COMMIT)
        if idx_ao:
            write_task.start()
        if idx_ai:
            read_task.start()

        # Starting the data acquisition/reproduction

        # Initialization

        if 'live' in args.plotting and idx_ai:live_Plot = myPlt.livePlot(args, number_of_channels_in, bufferSize, sample_rate, selection, copy.copy(chSelect))

        timeCounter = 0
        blockidx = 0
        values_read = np.zeros((sum(number_of_channels_in), int(number_of_samples + lastBufferZeroPadding)))
        current_buffer = np.zeros((sum(number_of_channels_in), int(bufferSize)))
        previous_buffer = np.zeros((sum(number_of_channels_in), int(bufferSize)))
        if sum(number_of_channels_in) == 1:
            current_buffer = np.array(current_buffer).reshape(1,len(current_buffer[0]))
            current = np.array(current_buffer).reshape(1,len(current_buffer[0]))
            previous_buffer = np.array(previous_buffer).reshape(1,len(previous_buffer[0]))

        # Main loop
        if idx_ai:
            pbar_ai = tqdm(total=number_of_samples)
            while timeCounter < number_of_samples:
                # The current read buffer
                current_buffer_raw = read_task.read(number_of_samples_per_channel=bufferSize)
                if sum(number_of_channels_in) > 1:
                    current_buffer = np.array(current_buffer_raw)
                    current = np.array(current_buffer_raw)
                else:
                    current_buffer[0,:] = np.array(current_buffer_raw)
                    current[0,:] = np.array(current_buffer_raw)

                if 'live' in args.plotting: live_Plot.livePlotClipp(current)
                # This is the variable that stores the data for saving
                values_read[:, timeCounter:timeCounter+bufferSize] = current_buffer
                # Calculations needed depending on the channel
                if signal_type[0]=='noise' and sum(number_of_channels_in) == 2 and selection != 'none':
                    previous_buffer, spectra, blockidx, H, _, HdB, H_phase, IR, gamma2, _ = TF.h1_estimator(current_buffer, previous_buffer, blockidx, calibrationData, micAmp, selection)
                    if 'live' in args.plotting: live_Plot.livePlotUpdate_H1(current, spectra, HdB, IR, gamma2)
                else:
                    if 'live' in args.plotting: live_Plot.livePlotUpdate(current)

                timeCounter += bufferSize
                pbar_ai.update(bufferSize)
            pbar_ai.close()
        elif idx_ao:
            time.sleep(sim_time)

    # Updating the dictionary with the measured data
    ch_count = 0
    for idx, device_idx in enumerate(idx_ai):
        if number_of_channels_in[idx] == 0 or idx_ai == []:
            continue
        for ch_num in range(number_of_channels_in[idx]):
            measurements.update({(channel_list[device_idx]['ai'][ch_num]): values_read[ch_count, :-lastBufferZeroPadding]})
            ch_count += 1

    if 'noise' in signal_type and number_of_channels_in[0] == 2 and selection != 'none':
        #Time and frequency vectors
        tVec = np.linspace(0, bufferSize / sample_rate, bufferSize)
        fftfreq = np.fft.rfftfreq(bufferSize, 1 / sample_rate)

        measurements.update({"Hij": H, "HdBij": HdB, "H_phase_ij": H_phase, "IRij": IR, "gamma2ij": gamma2, "fftfreq": fftfreq, "tVec": tVec})

    # Time signals plotting
    if 'timeSig' in args.plotting:
        time_signals_plot = {'xAxisTitle': 'time in s',
                             'yAxisTitle': 'Amplitude',
                             'plotTitle':   args.save_file + ' time_signals',
                             'scale':'lin'}

        full_tVec = np.linspace(0, sim_time, int(len(values_read[0,:])))
        time_signals = values_read
        if refChannel !='none':
            time_signals[:selection, :] = (values_read[:selection, :] / micAmp) * calibrationData[0]
            time_signals[selection+1:, :] = (values_read[selection+1:, :] / micAmp) * calibrationData[0]
        else:
            time_signals = (values_read / micAmp) * calibrationData[0]

        for idx in range(sum(number_of_channels_in)):
            time_signals_plot.update({str(idx):[full_tVec, time_signals[idx,:],"Channel:" + chSelect[idx]]})

        myPlt.singlePlot(time_signals_plot)


    return measurements
