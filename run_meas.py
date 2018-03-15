#Python libraries
import argparse
import glob
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import sys
import os, errno
import numpy as np
import time as ti
import scipy.io as sio
#Program libraries
import ni_tools as ni
import post_processing as pp
import meas_signal as sig
import utils as gz

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
    parser.add_argument("-pp", "--postProcess", help="Post proccessing type", type=str, choices={"no","TF"} )
    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 0.4')

    # Parse arguments
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parseArguments()
    #Setting the default values for the input arguments
    if not args.signalType: args.signalType = "pink_noise"
    if not args.pad_samples: args.pad_samples = 5000
    if not args.channelsIn: args.channelsIn = [0]
    if not args.channelsOut: args.channelsOut = [0]
    if not args.aiRange: args.aiRange = 5
    if not args.aoRange: args.aoRange = 1
    if not args.save_file: args.save_file = "no"
    if not args.bufferSize: args.bufferSize = 8192
    if not args.postProcess or args.calibration: args.postProcess = "no"
    if not args.micAmplification: args.micAmplification = 1
    if args.postProcess != "no": args.time == args.sampleRate == 0
    if not args.calibration: args.calibration = "no"

    #figure out the saved filename
    if args.save_file != "no": filename = args.save_file + "_"
    elif args.calibration != "no": filename = args.calibration 
    else: filename = "measurement_"

    #Configuring the setup for the calibration measurement
    if args.calibration =="new":
        args.channelsIn = [1]
        args.channelsOut = [0]
        args.pad_samples = 5000
        args.save_file = "no"
        args.bufferSize = 8192
    if args.postProcess != "no": args.save_file = filename

    #Printing the final arguements
    print("You are running the script with arguments: ")
    for a in args.__dict__:
        print(str(a) + ": " + str(args.__dict__[a]))
    #Creating the directories for saving the data if it does not exist
    directory = "aquired_data\\"
    meas_directory = "aquired_data\\measurement_"
    Caldirectory = "aquired_data\\calibration_files"
    currentDTime = ti.strftime('%y%m%d-%H%M')
    currentDay = ti.strftime('%y%m%d')
    gz.create_dir(directory)
    gz.create_dir(Caldirectory)

    #Calibration. New calibration measurement or loading the calibration data
    if args.calibration =="new":
        cal_filename = input("Please enter the name for the calibration file:")
        if not cal_filename: cal_filename = filename
        cal_postFilename = Caldirectory + "\\" + cal_filename + "_Cal"
        sensitivity = input("Please insert microphone sensitivity in mV/Pa(default 47.1 mV/Pa):")
        if not sensitivity: sensitivity = 47.1
        meas = ni.ni_io_tf(args)
        calibrationData = pp.mic_calibration(meas, sensitivity)
        np.save(cal_postFilename, calibrationData)
        sys.exit()
    elif args.calibration != "new" and args.calibration !="no":
        calFilename = Caldirectory + "\\" + args.calibration + "_Cal.npy"
        calibrationData = np.load(calFilename)
    else:
        print("No calibration data given.")
        calibrationData = [1, 1]


    #Post processing
    if args.postProcess != "no":
        #Choosing the directory
        print("\nPlease select a directory")
        idx = 0
        dir_names = []
        for name in os.listdir(path="aquired_data"):
            if os.path.isdir(directory + name):
                dir_names.append(name)
                print("[" + str(idx) + "]" + name)
                idx += 1
        selection = input()
        if not selection: selection = 0
        selected_directory = directory + dir_names[int(selection)]
        #Create the directories if not pressent
        TFdirectory = selected_directory + "\\TFs"
        gz.create_dir(TFdirectory)
        matTF_directory = selected_directory + "\\matLab\\TFs"
        gz.create_dir(matTF_directory)

        #File selection
        print("\nPlease enter the name of the file with the raw data: ")
        idx = 0
        file_names = []
        print(selected_directory)
        for name in sorted(glob.glob(selected_directory + '/*.np[yz]')):
                if "_Cal" in name: continue
                file_names.append(name)
                print("[" + str(idx) + "]" + name.rsplit('\\',1)[1].rsplit('.',1)[0])
                idx += 1
        selection = input()
        #Setting up the initial parameters
        blockSize = input("\nPlease enter desired blocksize for the analysis: ")
        if not blockSize: blockSize = 8192
        #If "all" is given processes all data in the specified directory
        if not selection:
            numpy_name = []
            numpy_vars = {}
            for filename in sorted(glob.glob(selected_directory + '/*.np[yz]')):
                if "_Cal" in filename: continue
                print("Processing file: " + filename)
                processedData = pp.h1_estimator(filename,  blockSize, calibrationData)
                #Saving the data
                postFilename = TFdirectory + "\\" + filename.rsplit('\\',1)[1].rsplit('.',1)[0] + "TFs_" + currentDTime
                matTF_postFilename = matTF_directory + "\\" + filename.rsplit('\\',1)[1].rsplit('.',1)[0] + "TFs_" + currentDTime
                np.save(postFilename, processedData)
                matsave = gz.saveToMat(processedData)
                sio.savemat(matTF_postFilename,{"object": matsave})
        else:
            filename = file_names[int(selection)]
            postFilename = filename
            processedData = pp.h1_estimator(postFilename,  blockSize, calibrationData)
            postFilename = TFdirectory + "\\" + filename.rsplit('\\',1)[1].rsplit('.',1)[0] + "TFs_" + currentDTime
            matTF_postFilename = matTF_directory + "\\" + filename.rsplit('\\',1)[1].rsplit('.',1)[0] + "TFs_" + currentDTime
            np.save(postFilename, processedData)
            matsave = gz.saveToMat(processedData)
            sio.savemat(matTF_postFilename,{"object": matsave})
    else:

    #New measurement
        meas = ni.ni_io_tf(args, calibrationData)
        # input("press enter")
        #Saving the new measurement
        if args.save_file!="no":
            meas_directory += currentDay
            gz.create_dir(meas_directory)
            postFilename = meas_directory + "\\" + args.save_file + "_" + currentDTime
            np.save(postFilename, meas)
            matsave = gz.saveToMat(meas)
            mat_directory = meas_directory + "\\matLab"
            gz.create_dir(mat_directory)
            mat_postFilename = mat_directory + "\\" + args.save_file + currentDTime
            sio.savemat(mat_postFilename,{"object": matsave})

    if sys.flags.interactive != 1 or not hasattr(QtCore, 'PYQT_VERSION'):
        pg.QtGui.QApplication.exec_()
