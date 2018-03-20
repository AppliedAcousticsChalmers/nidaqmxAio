# Python libraries
import argparse
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import sys
import numpy as np
# Program libraries
import systemUtils as sy
import ni_tools as ni
import post_processing as pp


def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("time", help="Measurement time in s", type=float)
    parser.add_argument("sampleRate", help="Measurement sample rate in Hz", type=float)

    # Optional arguments
    parser.add_argument("-tp", "--signalType", help="Type of test signal.", type=str, choices={"white_noise", "pink_noise", "debug"} )
    parser.add_argument("-pad", "--pad_samples", help="Number of zeros added to the end of the signal.", type=int)
    parser.add_argument("-ai", "--channelsIn", help="Number of channels per input module", nargs='+', type=int)
    parser.add_argument("-ao", "--channelsOut", help="Number of channels per output module", nargs='+', type=int)
    parser.add_argument("-aip", "--aiRange", help="Analog input absolute peak in Volts", type=float)
    parser.add_argument("-aop", "--aoRange", help="Analog output absolute peak in Volts", type=float)
    parser.add_argument("-sv", "--save_file", help="Exported data to filename.", type=str)
    parser.add_argument("-bf", "--bufferSize", help="Input buffer size in samples", type=int)
    parser.add_argument("-cal", "--calibration", help="typing \"new\" runs a new calibration measurement, otherwise specify the calibration file", type=str)
    parser.add_argument("-micA", "--micAmplification", help="Specify the amplification factior of he microphone", type=float)
    parser.add_argument("-pp", "--postProcess", help="Post proccessing type", type=str, choices={"no","TF","T60"} )
    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 0.4')

    # Parse arguments
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = parseArguments()
    # Setting the default values for the input arguments
    if not args.signalType: args.signalType = "pink_noise"
    if not args.pad_samples: args.pad_samples = 5000
    if not args.channelsIn: args.channelsIn = [0]
    if not args.channelsOut: args.channelsOut = [0]
    if not args.aiRange: args.aiRange = 5
    if not args.aoRange: args.aoRange = 1
    if not args.save_file: args.save_file = "no"
    if not args.bufferSize: args.bufferSize = 8192
    if not args.postProcess : args.postProcess = "no"
    if not args.micAmplification: args.micAmplification = 1
    if args.postProcess != "no": args.time == args.sampleRate == 0
    if not args.calibration: args.calibration = "no"

    # figure out the saved filename
    if args.save_file != "no": filename = args.save_file + "_"
    elif args.calibration != "no": filename = args.calibration
    else: filename = "measurement_"

    # Configuring the setup for the calibration measurement
    if args.calibration =="new":
        args.channelsIn = [1]
        args.channelsOut = [0]
        args.pad_samples = 5000
        args.save_file = "no"
        args.bufferSize = 8192
    if args.postProcess != "no": args.save_file = filename

    # Printing the final arguements
    print("You are running the script with arguments: ")
    for a in args.__dict__:
        print(str(a) + ": " + str(args.__dict__[a]))
    # Creating the directories for saving the data if it does not exist
    directory = "acquired_data"
    meas_directory = "acquired_data\\measurement_"
    Caldirectory = "acquired_data\\calibration_files"
    sy.create_dir(directory)
    sy.create_dir(Caldirectory)

    # Calibration. New calibration measurement or loading the calibration data
    if args.calibration == "new":
        cal_filename = input("Please enter the name for the calibration file:")
        if not cal_filename: cal_filename = filename
        cal_postFilename = Caldirectory + "\\" + cal_filename + "_Cal"
        sensitivity = input("Please insert microphone sensitivity in mV/Pa(default 47.1 mV/Pa):")
        if not sensitivity: sensitivity = 47.1
        meas = ni.ni_io_tf(args)
        calibrationData = pp.mic_calibration(meas, sensitivity)
        np.save(cal_postFilename, calibrationData)
        sys.exit()
    elif args.calibration != "new" and args.calibration != "no":
        calFilename = Caldirectory + "\\" + args.calibration + "_Cal.npy"
        calibrationData = np.load(calFilename)
    else:
        print("No calibration data given.")
        calibrationData = [1, 1]

    # Post processing
    if args.postProcess != "no":
        # Choosing the directory
        selected_diretory = sy.dir_select(directory)
        print(selected_diretory)
        # File selection
        filenames = sy.file_select(selected_diretory)
        # Setting up the initial parameters
        for current_file in filenames:
            print("Processing file: " + current_file)
            if args.postProcess == "TF": processedData = pp.h1_estimator(current_file,  args.bufferSize, calibrationData)
            # Saving the data
            sy.file_save(processedData, current_file, selected_diretory, args.postProcess)
    else:

        # New measurement
        meas = ni.ni_io_tf(args, calibrationData)
        # Saving the new measurement
        if args.save_file != "no":
            sy.file_save(meas, args.save_file, meas_directory)

    if sys.flags.interactive != 1 or not hasattr(QtCore, 'PYQT_VERSION'):
        pg.QtGui.QApplication.exec_()
