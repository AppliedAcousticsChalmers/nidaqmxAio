import os, errno
import glob
import time as ti
import scipy.io as sio
from collections.abc import Mapping
import numpy as np


# Prints the sub directories for given directory, asks the user to select one and returns the selected path name
def dir_select(directory):
    print("\nPlease select a directory")
    idx = 0
    dir_names = []
    for name in os.listdir(path=directory):
        if "calibration_files" in name:continue
        if os.path.isdir(directory + "\\" + name):
            dir_names.append(name)
            print("[" + str(idx) + "]" + name)
            idx += 1
    selection = input()
    if not selection: selection = 0
    selected_directory = directory + "\\" + dir_names[int(selection)] + "\\"

    return selected_directory


# Replaces the keys in a dictionary
def rec_key_replace(obj):
    if isinstance(obj, Mapping):
        return {key.replace('/', '_'): rec_key_replace(val) for key, val in obj.items()}
    return obj


# Takes a dictonary and formats it so it can be saved to matLab correctly
def saveToMat(obj):
    # obj = obj.item()
    mat = rec_key_replace(obj)
    return mat


# Creates a directory if it does not exist
def create_dir(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    return


# Prints the .np[zy] files pressent in the specified directory and returns the file name selected by the user.
# If no selection, return the first filename of the list.
# if 'all', return a list with all the filenames.
def file_select(directory):
        # File selection
        print("\nPlease enter the name of the file with the raw data: ")
        idx = 0
        file_names = []
        print(directory)
        for name in sorted(glob.glob(directory + '/*.np[yz]')):
                if "_Cal" in name: continue
                file_names.append(name)
                print("[" + str(idx) + "]" + name.rsplit('\\', 1)[1].rsplit('.', 1)[0])
                idx += 1
        selection = input()

        if not selection: return [file_names[0]]
        elif selection == "all": return file_names
        else: return [file_names[int(selection)]]


def file_save(data, filename, directory, extention=[]):
    currentDTime = ti.strftime('%y%m%d-%H%M')
    currentDay = ti.strftime('%y%m%d')
    # For the new measurement
    if not extention:
        # Main measurement directory creation
        directory += currentDay
        create_dir(directory)
        # Filename creation for the current meassurement
        postFilename = directory + "\\" + filename + "_" + currentDTime
        np.save(postFilename, data)
        # MatLab directory creation and export file
        matsave = saveToMat(data)
        mat_directory = directory + "\\matLab"
        create_dir(mat_directory)
        mat_filename = mat_directory + "\\" + filename + currentDTime
        sio.savemat(mat_filename, {"object": matsave})
    # For the postProcess scripts
    else:
        # Directory creation
        directory += extention + "\\"
        create_dir(directory)
        mat_directory = directory + "\\matLab_" + extention
        create_dir(mat_directory)
        filename = directory + "\\" + filename.rsplit('\\', 1)[1].rsplit('.', 1)[0] + "_" + extention + "_" + currentDTime
        mat_filename = mat_directory + "\\" + filename.rsplit('\\', 1)[1].rsplit('.', 1)[0] + "_" + extention + "_" + currentDTime
        np.save(filename, data)
        matsave = saveToMat(data)
        sio.savemat(mat_filename, {"object": matsave})

        return
