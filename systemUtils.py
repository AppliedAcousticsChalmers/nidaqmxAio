import errno
import os
import os.path
from pathlib import Path
from IPython.display import clear_output
import glob
import time as ti
import scipy.io as sio
from collections.abc import Mapping
import numpy as np
import math
import sys
from pyqtgraph.Qt import QtCore
import pyqtgraph as pg


import ni_tools as ni
import post_processing as pp


def fileSystem(directory, preSelect=''):
    '''
    Prints the sub directories for given directory and the .npy, .npz files,
    asks for user input.
    The names in the excluded_options list are excluded from the function

    User inputs:

    - If a number corresponding to a file is given, the function returns the
    path name of the file.
    - If a number corresponding to a directory is given, the function open the
    directory and asks again for
    user input.
    - If a number not corresponding to a file or a directory is given, the user
    is asked to give input again.
    - If all is given the function returns a list with all the .npy or .npz
    file paths.
    - If cd. is given the function goes up one folder and asks for user input.
    '''

    # excluded_options = ['calibration_files', 'matLab', 'matLab_TF',
    #                     'matlab_T60', 'TF', 'T60', "_Cal"]

    print("\nPlease select an option:")
    while True:
        try:
            idx = 0
            dir_names = []
            files_in_folder = \
                [fn for fn in glob.glob(directory + '/*.np[yz]')
                 if not os.path.basename(fn).endswith('_Cal.npy')]
            clear_output()
            for name in os.listdir(path=directory):
                # if name in excluded_options: continue
                if (directory + r"/" + name) in files_in_folder:
                    print("[" + str(idx) + "]" + name)
                elif os.path.isdir(directory + r'/' + name):
                    print("[" + str(idx) + "]" + name)
                else:
                    continue
                dir_names.append(name)
                idx += 1

            # Selecting the first option if no input is given
            if preSelect:
                selection = preSelect
            else:
                selection = input()
                if not selection:
                    selection = 0

            # Checking if the selection is a directory
            if selection != "cd." and selection != "all":
                is_directory = os.path.isdir(directory + r'/' +
                                             dir_names[int(selection)])

            if selection == "all":
                output = files_in_folder, directory
                break
            elif selection == 'cd.':
                directory = directory.rsplit(r"/", 1)[0]
                pass
            elif is_directory is True:
                directory = directory + r'/' + dir_names[int(selection)]
            else:
                selection = int(selection)
                output = [directory + r"/" + dir_names[selection]], directory
                break
            print("")
            print('Current directory:')
            print(directory)
            print("")

        except Exception:
            print('Please select one of the following')

    return output


def listToStr(data):
    output = []
    for i in data:
        string = str(i)
        output.append(string)
    return output


def rec_key_replace(obj):
    '''
    Replaces the given character in an object
    '''
    if isinstance(obj, Mapping):
        return {key.replace('/', '_'): rec_key_replace(val)
                for key, val in obj.items()}
    return obj


def saveToMat(obj):
    '''
    Formats the dictionary in a matLab compatible manner
    '''
    mat = rec_key_replace(obj)
    return mat


def create_dir(directory):
    '''
    Creates a directory if it does not exist
    '''
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    return


def file_save(data, filename, directory, extension=[], options=[]):
    '''
    data - data to be saved
    filename - the filename of the output filename
    directory - the directory to save the file into
    extension - creates a new folder inside the directory with the name given
    in the extension and adds the extension to the filename
    options - outputs the options .txt file
    '''
    currentDTime = ti.strftime('%y%m%d-%H%M')
    currentDay = ti.strftime('%y%m%d')

    # New measurement
    if not extension:
        # Main measurement directory creation
        directory += currentDay
        create_dir(directory)
        # Filename creation for the current meassurement
        postFilename = directory + r"/" + filename + "_" + currentDTime
        np.save(postFilename, data)
        # save configuration options
        f = open(postFilename + ".txt", "w")
        f.write(options)
        f.close()
        # MatLab directory creation and export file
        matsave = saveToMat(data)
        mat_directory = directory + r"/" + \
            directory.rsplit(r'/', 1)[1].rsplit('.', 1)[0] + "_" + "matLab"
        create_dir(mat_directory)
        mat_filename = mat_directory + r"/" + filename + "_" + currentDTime
        sio.savemat(mat_filename, {"object": matsave})
        f = open(mat_filename + ".txt", "w")
        f.write(options)
        f.close()

        return postFilename, directory

    # PostProcess scripts
    else:
        # Directory creation
        directory += r"/" + extension + r"/"
        create_dir(directory)
        mat_directory = directory + "matLab_" + extension
        create_dir(mat_directory)
        filename_out = directory + r"/" + \
            filename.rsplit('/', 1)[1].rsplit('.', 1)[0] + \
            "_" + extension + "_" + currentDTime
        mat_filename_out = mat_directory + r"/" + \
            filename.rsplit(r'/', 1)[1].rsplit('.', 1)[0] + \
            "_" + extension + "_" + currentDTime
        np.save(filename_out, data)
        matsave = saveToMat(data)
        sio.savemat(mat_filename_out, {"object": matsave})

        return


def simple_file_save(data, filename, directory):
    '''
    Same as the file_save function but without
    the extension and the options variables
    '''
    # Create directory if it does not exist
    create_dir(directory)
    # Create the filename using the input given so
    # it's understandable by the OS
    postFilename = directory + r"/" + filename
    # Save file in .npy format
    np.save(postFilename, data)
    # Create the matLab sub directory name
    mat_directory = directory + r"/matLab/"
    # Create matLab directory
    create_dir(mat_directory)
    # Create the filename using the input given so
    # it's understandable by the os
    mat_filename = mat_directory + filename
    # Format the data to be combatible with matLab
    matsave = saveToMat(data)
    # Save the .mat file
    sio.savemat(mat_filename, {"object": matsave})

    return


def sampleRateCorrection(input_sample_rate, fm=13.1072e6,
                         use_only_int_sr=True):
    '''
    The NI 9234, 9260 modules have an internal clock set at fm=13.1072MHz.
    This way they support sample rates that are given by the formula:
    fs = (fm / 256) / n where is an integer ranging from 1 to 31.
    The function gets the input sample rate and outputs one of the supported fs
    of the modules.
    If use_only_int_sr is set to True, only the sample rates that are integer
    numbers are selected.
    If an external clock is used then the fm parameter should be set to that
    clock's sample value.
    '''
    # fs = lambda n, fm: (fm / 256) / n
    def fs(n, fm):
        return (fm / 256) / n
    supported_sr = []
    for n in range(1, 32):
        if use_only_int_sr is True:
            if fs(n, fm) == math.floor(fs(n, fm)):
                supported_sr.append(fs(n, fm))
        else:
            supported_sr.append(fs(n, fm))
    supported_sr = np.array(supported_sr)
    selection = np.argmin(abs(supported_sr - input_sample_rate))
    corrected_sample_rate = supported_sr[selection]

    return corrected_sample_rate


def calibrationFilesCheck(args, directory):
    '''
    Handles the calibration options for the nidaqmxAio through
    the command line.
    Inputs:
    - args: the parser arguments of the program
    - directory: the directory that the function will check for
    calibration files
    '''
    if args.calibration == 'uncalibrated':
        print('Calibration canceled by user.')
        calibrationData = [1, 1]
        return calibrationData, args

    if args.calibration:
        cal_postFilename = directory + "\\" + args.calibration + "_Cal"
        cal_path = Path(cal_postFilename + ".npy")
        if os.path.isfile(cal_path):
            calibrationData = np.load(cal_postFilename + '.npy',
                                      allow_pickle=True)
        else:
            print("Calibration file '" + str(cal_path) + "' doesn't exist.\n")
            user_input = input("Press < Enter > to continue without " +
                               "calibration, or < c > to calibrate and save: ")
            if user_input == 'c':
                args.postProcess = []
                args.newMeasurement = 0
                meas = ni.ni_io_tf(args, cal=True)
                calibrationData = pp.mic_calibration(meas, args.sensitivity)

                simple_file_save(calibrationData, args.calibration + "_Cal",
                                 directory)
            else:
                print("No calibration data given")
                calibrationData = [1, 1]
    else:
        pressent_calibration_files = [fn for fn in glob.glob(directory +
                                                             '/*.np[yz]')
                                      if os.path.basename(fn)
                                      .endswith('_Cal.npy')]
        if pressent_calibration_files:
            print("\nCalibration files present in current directory:")
            while True:
                try:
                    idx = 0
                    dir_names = []
                    for name in os.listdir(path=directory):
                        if (directory + r"/" + name) \
                           in pressent_calibration_files:
                            print("[" + str(idx) + "]" + name)
                        else:
                            continue
                        dir_names.append(name)
                        idx += 1

                    print("[" + str(idx) + "]No calibration.")
                    selection = input()
                    if not selection:
                        selection = 0
                    selection = int(selection)
                    if selection == idx:
                        print("No calibration data given")
                        calibrationData = [1, 1]
                        return calibrationData, args

                    selected_calibration_file = directory + r"/" + \
                        dir_names[selection]
                    calibrationData = np.load(selected_calibration_file,
                                              allow_pickle=True)
                    break
                except Exception:
                    print('Please select one of the following')

        else:
            print("No calibration data given")
            calibrationData = [1, 1]

    return calibrationData, args
