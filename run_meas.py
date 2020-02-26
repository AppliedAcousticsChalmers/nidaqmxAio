# Python libraries
from pyqtgraph.Qt import QtCore
import pyqtgraph as pg
import sys
import configargparse

# Program libraries
import systemUtils as sy
import ni_tools as ni
import post_processing as pp


if __name__ == '__main__':
    p = configargparse.ArgParser(default_config_files=['default_config.yml'])
    p.add(
        '-c', '--config', is_config_file=True,
        help='config file path')
    p.add(
        "--time", required=True,
        help="Measurement time in s",
        type=int)
    p.add(
        "--sampleRate", required=True,
        help="Measurement sample rate in Hz",
        type=int)
    p.add(
        "--newMeasurement", required=True,
        help="Perform measurement",
        type=int, default=0, choices={0, 1})
    # Optional arguments
    p.add(
        "-sig", "--signalType",
        help="Type of test signal",
        nargs='+', default=['noise_pink'])
    p.add(
        "-pad", "--pad_samples",
        help="Number of zeros added to the end of the signal",
        type=int, default=512)
    p.add(
        "-ai", "--channelsIn",
        help="Number of channels per input module",
        nargs='+', type=int, default=[0])
    p.add(
        "-ao", "--channelsOut",
        help="Number of channels per output module",
        nargs='+', type=int, default=[0])
    p.add(
        "-aip", "--aiRange",
        help="Analog input absolute peak in Volts",
        type=float, default=5)
    p.add(
        "-aop", "--aoRange",
        help="Analog output absolute peak in Volts",
        type=float, default=0.001)
    p.add(
        "-sv", "--save_file",
        help="Exported data to filename.",
        type=str, default='measurement')
    p.add(
        "-dir", "--root_directory",
        help="Root directory, for processing and storing data.",
        type=str, default='acquired_data')
    p.add(
        "-bf", "--bufferSize",
        help="Input buffer size in samples",
        type=int, default=8192)
    p.add(
        "-cal", "--calibration",
        help="Specify the calibration file, if it does not exist,"
        + " ask if a new calibration measurement is needed",
        type=str)
    p.add(
        "-micA", "--micAmp",
        help="Specify the amplification factor of he microphone",
        type=float, default=1)
    p.add(
        "-sens", "--sensitivity",
        help="microphone sensitivity in mV/Pa",
        type=float, default=47.1)
    p.add(
        "-pp", "--postProcess",
        help="Post processing type. Do not define, if no post processing script is needed.",
        type=str, choices={"TF", "RAC"})
    p.add(
        "-cT", "--cutoffTime",
        help="Measurement time after the end of the signal, in s. Must be int",
        type=int, default=0)
    p.add(
        "-plt", "--plotting",
        help="Plots to display. Options: live, TF, timeSig, "
        + "T60_one_band, T60_3rd",
        nargs='+', default=['live'])
    p.add(
        "-fRange", "--frequencyRange",
        help="Frequency range for T60 calculations,"
        + " example [fmin, fmax, bandwidth]",
        nargs='+', default=[20, 10000, 'third'])
    p.add(
        "-refCh", "--refferenceChannel",
        help="Pre-specify which channel will be used as reference",
        default="")
    p.add(
        "-nt", "--note",
        help="Adds a text note to the save file.",
        type=str, default="")

    # Parse arguments
    args = p.parse_args()
    # Print used arguments
    print(p.format_values())
    # Create directories for saving data
    directory = args.root_directory
    meas_directory = directory + r"/measurement_"
    sy.create_dir(directory)

    # Calibration
    calibrationData, args = sy.calibrationFilesCheck(args, meas_directory)
    print("Calibration finished.")

    # New measurement
    if args.newMeasurement == 1:
        # Measurement
        meas = ni.ni_io_tf(args, calibrationData)
        # Save the new measurement
        filenames, meas_directory = \
            sy.file_save(
                meas,
                args.save_file,
                meas_directory,
                options=p.format_values())

    # Post processing
    if args.postProcess:
        # File selection
        if args.newMeasurement == 0:
            # Choosing the directory and the files
            filenames, selected_directory = sy.fileSystem(directory)
            # Calibration
            calibrationData, args = \
                sy.calibrationFilesCheck(args, selected_directory)
        else:
            selected_directory = meas_directory
            filenames = [filenames + ".npy"]
        # Setting up the initial parameters
        for current_file in filenames:
            print("Processing file: " + current_file)
            if args.postProcess == "TF":
                processedData = pp.TFcalc(
                    current_file,
                    args.bufferSize,
                    calibrationData,
                    args.plotting)
            elif args.postProcess == "RAC":
                processedData = pp.RAC(
                    current_file,
                    args.bufferSize,
                    args.frequencyRange[0],
                    args.frequencyRange[1],
                    args.frequencyRange[2],
                    calibrationData,
                    args.plotting)
        # Saving the data
            sy.file_save(
                processedData,
                current_file,
                selected_directory,
                args.postProcess)

    if (sys.flags.interactive != 1 and not args.refferenceChannel) \
       or not hasattr(QtCore, 'PYQT_VERSION') \
       and 'live' in args.plotting:
        pg.QtGui.QApplication.exec_()
