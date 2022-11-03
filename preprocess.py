''' 
Last changes made on: 2022-Nov-02
@author: Mateusz Fido, mateusz.fido@org.chem.ethz.ch
ETH Zürich

This script reads in and parses .mzml files using the Pyteomics
metabolomic library for Python. It calculates average mass spectra
by linearly interpolating their intensity on a resampled m/z axis 
with given resolution, then writes them as .csv files for purposes 
of graphing and/or analysis with other software. 

It also performs pick peaking by using scipy function find_peaks using
height cutoff criteria and stores peak data that can then be passed
onto functions calculating time traces (XICs) of these peaks.

TODO: correlate the peak time traces with the TIC ? 

References:
===========
Matlab Code of Jiayi Liang and Miguel de Figueiredo 
Python Code of Cedric Wüthrich
https://pyteomics.readthedocs.io/en/latest/index.html

'''

from math import ceil
import os, time
from pathlib import Path
import numpy as np
from pyteomics import mzml
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt

st = time.time()    # log execution time

'''Hardcoded data paths for debugging'''
#
# Windows data path: 
# PATH = Path('C:\\Users\\drago\\OneDrive - ETH Zurich\\Kispi\\data-analysis\\test_mzml')
#
# MacOS data path: 
PATH = Path('/Users/mateuszfido/OneDrive - ETH Zurich/Kispi/data-analysis/test_mzml/')


'''Constants used for the new mz_axis'''
MZ_MIN = 50
MZ_MAX = 500
RES = 0.005
DATA_POINTS = int((MZ_MAX - MZ_MIN)/RES)

MZ_AXIS = np.linspace(MZ_MIN, MZ_MAX, DATA_POINTS, dtype=float)

def read_mzml(path = PATH):
    '''
    Parses .mzml data specified in the input path. \n
    -----------------------------------------------------------------
    Input: Path to a directory containing .mzml files. \n
    -----------------------------------------------------------------
    Returns: Python List of File objects.
    '''

    folder_contents = os.listdir(path)

    filelist = []   # Create empty placeholder list 

    for file in folder_contents:
        if file.endswith(".mzml"):
            filelist.append(file)   # Append filenames

    return filelist

def read_csv(path = PATH):
    '''
    Parses .csv data specified in the input path. \n
    -----------------------------------------------------------------
    Input: Path to a directory containing .csv files. \n
    -----------------------------------------------------------------
    Returns: Python List of File objects.
    '''
    path_csv = path / "Extracted_Data/"
    folder_contents = os.listdir(path_csv)

    filelist_csv = []   # Create empty placeholder list 

    for file in folder_contents:
        if file.endswith(".csv"):
            filelist_csv.append(file)   # Append filenames

    return filelist_csv

def average():
    '''
    Reads .mzml files using read_mzml. Performs averaging of the mass spectra (scans)
    by linearly interpolating their intensities over a new, resampled, linearly spaced
    m/z axis with resolution specified as a constant at the top of the script.
    Saves the averaged mass spectra to .csv files in two comma-separated columns: 
    m/z and intensity.
    Memory-loading at high resolution.
    ------
    Input: None \n
    ------
    Returns: None \n'''

    filelist = read_mzml()

    path_to_save = PATH / "Extracted_Data/"

    try:
        os.mkdir(path_to_save)
    except:
        print("Directory already exists")   # Check if folder already exists, if not, create it

    for i in range(0, len(filelist)):   # Iterate over all files in path
        curr_file = PATH / filelist[i]
        mzml_file = mzml.MzML(str(curr_file))  # Instantiate the MzML reader object

        intensities = np.zeros(len(MZ_AXIS), dtype=float)
        background = np.zeros(len(MZ_AXIS), dtype=float)

        path_to_save_parsed = path_to_save / "{}avg.csv".format(filelist[i]).replace('.mzml', "_")

        with open(path_to_save_parsed, "w+") as average_csv:
            for j in range(0, len(mzml_file)):
                mz_array = np.ndarray.tolist(mzml_file[j]['m/z array'])
                intensity_array = np.ndarray.tolist(mzml_file[j]['intensity array'])
                int_interp = np.interp(MZ_AXIS, mz_array, intensity_array, left = 0, right = 0)
                intensities += int_interp
                
                while j < 20:   # Consider first 20 scans as background 
                    background += int_interp
                    break

            avg_intensity = intensities / len(mzml_file)   # Average the signal
            avg_background = background / 20                # Average the background 
            corr_intensity = avg_intensity - avg_background # Subtract averaged background signal
           
            data_matrix = np.array((MZ_AXIS, corr_intensity), dtype=float)

            data_matrix = np.transpose(data_matrix)                                 # Save background-corrected, 
            np.savetxt(average_csv, data_matrix, delimiter=",", encoding='utf-8')   # resampled m/z-intensities


class Peak():
    '''Support class for Peak objects. \n
    -----------
    Parameters:
    -----------
    mz: array of mz values computed from indices of intensity axis, returned by scipy.find_peaks()
    ints: integer, intensity value of the picked peak
    width: non-zero integer, peak width calculated by scipy.find_peaks()'''

    def __init__(self, mz, ints, width):
        self.mz = mz
        self.ints = ints
        self.width = width

    def __str__(self):
        return f"{self.mz}: {self.ints}"

def peak_pick():
    '''Performs peak-picking and integration of peaks from average mass spectra
    by using scipy.find_peaks() (finds local maxima by comparison of nearest points).\n
    
    Input: None\n
    
    Output: List of Peak objects'''

    avg_filelist = read_csv()

    peaklist = []

    for file_csv in avg_filelist:
        with open(PATH / "Extracted_Data/" / file_csv, "r") as curr_file:
            data = np.genfromtxt(curr_file, delimiter=',')  # Convert .csv into a numpy ndarray object
            corr_intensity = data.transpose()[1]            # Transpose it into rows and get the intensity array
            peaks = find_peaks(corr_intensity, height = 100)    # Find peaks by cutting off at a given intensity to remove noise
            
            # Find full widths at half maximum (FWHM) for later integration
            widths, width_heights, left_fwhm, right_fwhm = peak_widths(corr_intensity, peaks[0], rel_height=0.5) 

            for l in range(0, len(peaks[0])):   # For all peaks found, extract their properties and append to List
                mz = MZ_AXIS[peaks[0][l]]
                ints = peaks[0][l]
                width = ceil(widths[l]/2)
                peak = Peak(mz, ints, width)
                peaklist.append(peak)

    return peaklist            

    '''Graphing interface for debugging purposes'''

    # print(peaks)
    # fig = plt.figure()
    # ax = fig.subplots()
    # plt.plot(mz_axis, corr_intensity)
    # ax.scatter(mz_axis[peaks[0]], peaks[1]['peak_heights'], color = 'r', s=15)
    # ax.legend()
    # ax.grid()
    # plt.show()

def time_trace(file):
    '''For each scan (i.e. time unit), iterates over every peak in the peaklist provided.
    Because of O(n^2) (or higher?) time complexity the function is quite resource-heavy. \n
    Input: .mzml file to iterate over \n
    Output: data: Numpy array of peaks x scans, where the first row comprises m/z value of each peak's apex. 
    tic: array of the total ion current in the measurement for further peak correlation.'''

    peaklist = peak_pick()

    scans = mzml.MzML(str(PATH / file))

    data = np.empty((len(scans)+1, len(peaklist)))   # Placeholder array for peaks x scans
    tic = np.empty(len(scans))                       # Placeholder array for the tic

    i = 1   # Create iterators for indexing of the output array
    j = 0

    for scan in scans:
        tic[i-1] = scan['total ion current']
        mz_array = np.ndarray.tolist(scan['m/z array'])
        intensity_array = np.ndarray.tolist(scan['intensity array'])
        int_interp = np.interp(MZ_AXIS, mz_array, intensity_array, left = 0, right = 0)
        data[i] = scan['index']
        for peak in peaklist:
            peak.intensities = int_interp[(peak.ints - peak.width) : (peak.ints + peak.width)]
            try:
                integral = np.trapz(peak.intensities,  dx = RES)
                data[i][j] = integral
            except(ValueError):
                data[i][j] = 0
                continue
            finally:
                data[0][j] = peak.mz
                j += 1
        i += 1
        j = 0

    return data, tic


def tic_correlation(data, tic):
    '''TODO'''
    filelist = read_mzml()
    counter = 0
    for file in filelist:
        counter += 1
        print(f"Processing file: {file}, {counter} out of {len(filelist)}")
        data, tic = time_trace(file)
    print("Finished processing timetraces.")

et = time.time()
elapsed_time = et - st
print("Execution time: ", round(elapsed_time, 2), " seconds.")  # Log script running time