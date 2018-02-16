import matplotlib.pyplot as plt
import nidaqmx.system
import collections
import numpy as np
import time
import colorednoise as cn
from nidaqmx.constants import AcquisitionType, TaskMode
from nidaqmx.stream_writers import AnalogMultiChannelWriter
from nidaqmx.stream_readers import AnalogMultiChannelReader

system = nidaqmx.system.System.local()


def create_signal(sample_rate, time_of_signal, pad_samples, signal_type, voltage_out):
    number_of_samples = int(time_of_signal*sample_rate)
    if signal_type == "pink_noise":
        signal = cn.powerlaw_psd_gaussian(1,number_of_samples)
    elif signal_type == "white_noise":
        signal = cn.powerlaw_psd_gaussian(0,number_of_samples)
    elif signal_type == "debug":
        signal = np.sin(np.linspace(0,89*np.pi,number_of_samples))*np.sin(np.linspace(0,2000*np.pi,number_of_samples))
    signal /= (np.max(abs(signal)))/voltage_out
    signal = np.pad(signal, (pad_samples, 0), 'constant', constant_values = [0, 0])

    return signal

def ni_io_tf(args):
    sim_time = args.time
    sample_rate = args.sampleRate
    signal = create_signal(args.sampleRate, args.time, args.pad_samples, args.signalType, args.aoRange)
    number_of_channels_in = args.channelsIn
    number_of_channels_out = args.channelsOut
    ai_range = args.aiRange
    ao_range = args.aoRange
    bufferSize = args.bufferSize

    number_of_samples = sim_time*sample_rate + args.pad_samples
    if not np.sum(number_of_channels_out)<= 1: signal = np.tile(signal, [np.sum(number_of_channels_out), 1])
    measurements = {'simulationTime': sim_time, 'SampleRate': sample_rate, 'Signal': args.signalType, 'Input range in Vrms': ai_range, 'Output range in Vrms': ao_range}

    channel_list = []
    for device in system.devices:
        channel_list.append({"ai":[channels_ai.name for _, channels_ai in enumerate(device.ai_physical_chans)],
                            "ao":[channels_ao.name for _, channels_ao in enumerate(device.ao_physical_chans)]})
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

        timeCounter = 0
        if idx_ai:
            while timeCounter <= number_of_samples:
                current_buffer = read_task.read(number_of_samples_per_channel=bufferSize)
                # print(np.shape(current_buffer),np.shape(values_read))
                values_read = np.append(values_read, current_buffer, axis=1*(sum(number_of_channels_in)!=1))
                timeCounter += bufferSize
        elif idx_ao:
            time.sleep(sim_time)

    ch_count = 0
    for idx, device_idx in enumerate(idx_ai):
        if number_of_channels_in[idx] == 0 or idx_ai==[] : continue
        for ch_num in range(number_of_channels_in[idx]):
            measurements.update({(channel_list[device_idx]['ai'][ch_num]):np.array(values_read)[ch_count, ...]})
            ch_count += 1
    # plt.plot(np.array(values_read))
    # plt.plot(signal[0, ...])
    plt.show()
    print(measurements)
    return measurements
