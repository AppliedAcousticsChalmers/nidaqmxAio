# Python libraries
import os.path
import argparse
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import sys
import numpy as np
from pathlib import Path
import configargparse

# Program libraries
import systemUtils as sy
import ni_tools as ni
import post_processing as pp


if __name__ == '__main__':
    p = configargparse.ArgParser(default_config_files=['default_config.yml'])
    p.add('-c', '--config', is_config_file=True, help='config file path')
    p.add("--time", required=True, help="Measurement time in s", type=int)
    p.add("--sampleRate", required=True, help="Measurement sample rate in Hz", type=int)
    p.add("--newMeasurement", required=True, help="Perform measurement", type=int, default=0, choices={0,1})
    # Optional arguments
    p.add("-sig", "--signalType", help="Type of test signal", nargs='+', default=['noise_pink'])
    p.add("-pad", "--pad_samples", help="Number of zeros added to the end of the signal",  type=int, default=512)
    p.add("-ai", "--channelsIn", help="Number of channels per input module", nargs='+', type=int, default=[0])
    p.add("-ao", "--channelsOut", help="Number of channels per output module", nargs='+', type=int, default=[0])
    p.add("-aip", "--aiRange", help="Analog input absolute peak in Volts", type=float, default=5)
    p.add("-aop", "--aoRange", help="Analog output absolute peak in Volts", type=float, default=0.001)
    p.add("-sv", "--save_file", help="Exported data to filename.", type=str, default='measurement')
    p.add("-bf", "--bufferSize", help="Input buffer size in samples", type=int, default=8192)
    p.add("-cal", "--calibration", help="Specify the calibration file, if it does not exist, ask if a new calibration measurement is needed", type=str)
    p.add("-micA", "--micAmp", help="Specify the amplification factior of he microphone", type=float, default=1)
    p.add("-sens", "--sensitivity", help="microphone sensitivity in mV/Pa", type=float, default=47.1)
    p.add("-pp", "--postProcess", help="Post processing type", type=str, choices={"TF", "RAC"})
    p.add("-cT", "--cutoffTime", help="Measurement time after the end of the signal, in s", type=int, default=0)
    p.add("-plt", "--plotting", help="Plots to display. Options: live, TF, timeSig, T60_one_band, T60_3rd", nargs='+', default=['live'])
    p.add("-fRange", "--frequencyRange", help="Frequency range for postProcess calculation example \[fmin, fmax, bandwidth\]", nargs='+', default=[20, 10000, 'third'])
    p.add("-refCh", "--refferenceChannel", help="Prespecify which channel will be used as reference", default = "")
    p.add("-cmt", "--comment", help="Adds a text comment to the save file.", type=str, default="")

    # Parse arguments
    args = p.parse_args()
    # Print used arguements
    print(p.format_values())
    # Create directories for saving data
    directory = "acquired_data"
    meas_directory = "acquired_data\\measurement_"
    Caldirectory = "acquired_data\\calibration_files"
    sy.create_dir(directory)
    sy.create_dir(Caldirectory)

    # Calibration. New calibration measurement or loading the calibration data
    if args.calibration:
        cal_postFilename = Caldirectory + "\\" + args.calibration + "_Cal"
        cal_path = Path(cal_postFilename + ".npy")
        if os.path.isfile(cal_path):
            calibrationData = np.load(cal_postFilename + '.npy')
        else:
            user_input = input("Calibration file '" + str(cal_path) + "' doesn't exist.\nPress <Enter> to continue without calibration, or <c> to calibrate and save:")
            if user_input == 'c':
                args.postProcess = []
                args.newMeasurement = 0
                meas = ni.ni_io_tf(args,cal=True)
                calibrationData = pp.mic_calibration(meas, args.sensitivity)

                sy.simple_file_save(calibrationData, args.calibration + "_Cal", Caldirectory)


            else:
                print("No calibration data given")
                calibrationData = [1, 1]
    else:
        print("No calibration data given")
        calibrationData = [1, 1]

    # New measurement
    if args.newMeasurement == 1:
        meas = ni.ni_io_tf(args, calibrationData)
        # Save the new measurement
        filenames, meas_directory = sy.file_save(meas, args.save_file, meas_directory, options=p.format_values())

    # Post processing
    if args.postProcess:
        # File selection
        if args.newMeasurement == 0:
            # Choosing the directory and the files
            filenames, selected_directory = sy.fileSystem(directory)
        else:
            selected_directory = meas_directory
            filenames = [filenames + ".npy"]
        # Setting up the initial parameters
        for current_file in filenames:
            print("Processing file: " + current_file)
            if args.postProcess == "TF":
                processedData = pp.TFcalc(current_file,  args.bufferSize, calibrationData, args.plotting)
            elif args.postProcess == "RAC":
                processedData = pp.RAC(current_file, args.bufferSize, args.frequencyRange[0], args.frequencyRange[1], args.frequencyRange[2], calibrationData, args.plotting)
        # Saving the data
            sy.file_save(processedData, current_file, selected_directory, args.postProcess)

    if (sys.flags.interactive != 1 and not args.refferenceChannel) or not hasattr(QtCore, 'PYQT_VERSION') and 'live' in args.plotting :
        pg.QtGui.QApplication.exec_()
