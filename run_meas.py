# Python libraries
import argparse
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import sys
import numpy as np
from pathlib import Path
# Program libraries
import systemUtils as sy
import ni_tools as ni
import post_processing as pp
import configargparse


if __name__ == '__main__':
    p = configargparse.ArgParser()
    p.add('-c', '--config', required=True, is_config_file=True, help='config file path')
    p.add("--time", required=True, help="Measurement time in s")
    p.add("--sampleRate", required=True, help="Measurement sample rate in Hz")
    # Optional arguments
    p.add("-sig", "--signalType", help="Type of test signal", nargs='+', default=['noise_pink'])
    p.add("-pad", "--pad_samples", help="Number of zeros added to the end of the signal", type=int, default=5000)
    p.add("-ai", "--channelsIn", help="Number of channels per input module", nargs='+', type=int, default=[0])
    p.add("-ao", "--channelsOut", help="Number of channels per output module", nargs='+', type=int, default=[0])
    p.add("-aip", "--aiRange", help="Analog input absolute peak in Volts", type=float, default=5)
    p.add("-aop", "--aoRange", help="Analog output absolute peak in Volts", type=float, default=[1])
    p.add("-sv", "--save_file", help="Exported data to filename.", type=str, default='measurement')
    p.add("-bf", "--bufferSize", help="Input buffer size in samples", type=int, default=8192)
    p.add("-cal", "--calibration", help="typing \"new\" runs a new calibration measurement, otherwise specify the calibration file", type=str)
    p.add("-micA", "--micAmp", help="Specify the amplification factior of he microphone", type=float, default=1)
    p.add("-sens", "--sensitivity", help="microphone sensitivity in mV/Pa", type=float, default=47.1)
    p.add("-pp", "--postProcess", help="Post proccessing type", type=str, choices={"TF", "T60"})
    p.add("-meas", "--measurement", required=True, help="Perform measurement", type=bool, default=True)
    # Print used arguements
    print(p.format_values)

    # Parse arguments
    args = p.parse_args()

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
        if cal_path.isfile():
            calibrationData = np.load(cal_postFilename)
        else:
            user_input = input("Calibration file '" + cal_path + "' doesn't exist.\nPress <Enter> to continue without calibration, or <c> to calibrate and save")
            if user_input == 'c':
                meas = ni.ni_io_tf(args)
                calibrationData = pp.mic_calibration(meas, args.sensitivity)
                np.save(cal_postFilename, calibrationData)
    else:
        print("No calibration data given")
        calibrationData = [1, 1]
    # Measurement
    if args.measurement == True:
        # New measurement
        meas = ni.ni_io_tf(args, calibrationData)
        # Save the new measurement
        sy.file_save(meas, args.save_file, meas_directory, options=p.format_values())

    # Post processing
    if args.postProcess:
        # Choosing the directory
        selected_diretory = sy.dir_select(directory)
        print(selected_diretory)
        # File selection
        if args.measurement == True:
            # FIX
            # post process only the data just measured
            # filenames = q
        else:
            filenames = sy.file_select(selected_diretory)

        # Setting up the initial parameters
        for current_file in filenames:
            print("Processing file: " + current_file)
            if args.postProcess == "TF":
                processedData = pp.h1_estimator(current_file,  args.bufferSize, calibrationData)
            # Saving the data
            sy.file_save(processedData, current_file, selected_diretory, args.postProcess)

    if sys.flags.interactive != 1 or not hasattr(QtCore, 'PYQT_VERSION'):
        pg.QtGui.QApplication.exec_()
