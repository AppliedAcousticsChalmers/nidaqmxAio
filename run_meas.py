#Python libraries
import argparse
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
    parser.add_argument("-sv", "--save_file", help="Option for exporting data to file.", type=bool)
    parser.add_argument("-bf", "--bufferSize", help="Input buffer size in samples", type=int)
    parser.add_argument("-cal", "--calibration", help="Run calibration measurement", type=bool)
    parser.add_argument("-pp", "--postProcess", help="Post proccessing type", type=str, choices={"no","TF"} )
    parser.add_argument("-opp", "--onlyPostProcess", help="Run the selected post processing, without measurement ", type=str, choices={"no","TF"} )
    parser.add_argument("-cd", "--calData", help="calibration data", type=float)
    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 0.0')

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
    if not args.save_file: args.save_file = False
    if not args.bufferSize: args.bufferSize = 8192
    if not args.postProcess or args.calibration==True: args.postProcess = "no"
    if not args.onlyPostProcess or args.calibration==True: args.onlyPostProcess = "no"
    if args.onlyPostProcess != "no": args.time == args.sampleRate == 0
    if not args.calibration: args.calibration = False
    if not args.calData: args.calData = [1, 1]

    #Configuring the setup for the calibration measurement
    if args.calibration ==True:
        args.channelsIn = [1]
        args.channelsOut = [0]
        args.pad_samples = 5000
        args.save_file = False
        args.bufferSize = 8192
    if args.postProcess != "no": args.save_file = True


    #Ensuring that only one type of postprocessing is taking place
    if args.onlyPostProcess != "no" and args.postProcess != "no":
        print("Please select whether a new measurement should be taken [y/n]")
        answer = input()
        if answer == "y" or "Y" or "Yes" or "yes":
            args.onlyPostProcess == "no"
        elif answer == "n" or "N" or "No" or "no":
            args.postProcess == "no"
            args.time == args.sampleRate == 0

    #Printing the final arguements
    print("You are running the script with arguments: ")
    for a in args.__dict__:
        print(str(a) + ": " + str(args.__dict__[a]))
    directory ="aquired_data\\"
    currentDTime = ti.strftime('%y%m%d-%H%M')

    #The calibration measurement
    if args.calibration==True:
        filename = input("Please enter the desired name, for the calibration file: ")
        if not filename: filename = "calibration"
        cal_postFilename =  filename + "_Cal"
        sensitivity = input("Please insert microphone sensitivity in mV/Pa(default 47.1 mV/Pa):")
        if not sensitivity: sensitivity = 47.1
        meas = ni.ni_io_tf(args)
        calibrationData = pp.mic_calibration(meas, sensitivity)
        np.save(cal_postFilename, calibrationData)

    #Post processing without new measurement
    elif args.onlyPostProcess != "no":
        calFilename =  input("Please enter the name of the file with the calibration data: ")
        post_calFilename =  calFilename + "_Cal"
        if not calFilename:
            print("No calibration data given.")
            calibrationData = [1, 47.1]
        else:
            calibrationData = np.load(post_calFilename + ".npy")

        rawDataFilename =  input("Please enter the name of the file with the raw data: ")
        postRawDataFilename = directory + rawDataFilename
        blockSize = input("Please enter desired blocksize for the analysis: ")
        if not blockSize: blockSize = 8192
        processedData = pp.h1_estimator(postRawDataFilename,  blockSize, calibrationData)
        postFilename = postRawDataFilename + "_TFs_" + currentDTime
        np.save(postFilename, processedData)
        matsave = gz.saveToMat(processedData)
        sio.savemat(postFilename,{"object": matsave})
    else:

    #New measurement
        meas = ni.ni_io_tf(args)
        #Saving the new measurement
        if args.save_file==True:
            filename =  input("Please enter the desired name, for the output file: ")
            if not filename: filename = 'meas_'
            postFilename = directory + filename + currentDTime
            np.save(postFilename, meas)
            matsave = gz.saveToMat(meas)
            sio.savemat(postFilename,{"object": matsave})

    #Post processing using the new measurement
    if args.postProcess == "TF":
        calFilename =  input("Please enter the name of the file with the calibration data: ")
        post_calFilename =  calFilename + "_Cal"
        if not calFilename:
            print("No calibration data given (using default sensitivity 47.1mV/Pa).")
            calibrationData = [1, 47.1]
        else:
            calibrationData = np.load(post_calFilename + ".npy")

        blockSize = input("Please enter desired blocksize for the analysis: ")
        if not blockSize: blockSize = 8192
        processedData = pp.h1_estimator(postFilename,  blockSize, calibrationData)
        postFilename = directory + filename + "_TFs_" + currentDTime
        np.save(postFilename, processedData)
        matsave = gz.saveToMat(processedData)
        sio.savemat(postFilename,{"object": matsave})
