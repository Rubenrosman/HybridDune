# ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    DEFINE INPUT/OUTPUT FILES   ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    
input_folder = r'O:\HybridDune experiment\data RBR, OSSI\copy RBR Udrive series1\QC'
input_file = 'S4P3 RBR2 p_rel - period 1.nc'
output_folder = r'O:\HybridDune experiment\data RBR, OSSI\copy RBR Udrive series1\processed'
output_file = 'Pressure sensor S4P3 RBR2 processed data - period 1.nc'

# IMPORT MODULES ===================================================================================================================
import os
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import xarray as xr
import matplotlib.pyplot as plt
import sys
from scipy.signal import welch
sys.path.append(r'C:\Users\dpoppema\Documents\GitHub\HybridDune\Ruben\Pressure_sensors\S1\RBR_05') # to find the puv.py file
import puv 

# ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    DEFINE FUNCTION FOR WAVE PROCESSING   ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
def wave_statistics_netcdf(input_folder, input_file, output_folder, output_file):
    """
    Function to process wave data from a netCDF file, resample it to 2D blocks, and calculate wave statistics.
    The function reads the input data, reshapes it, calculates wave statistics using puv.py and saves the processed data to a new netCDF file.
    """
    # INPUT/OUTPUT FILES, SETTINGS ======================================================================================================
    # Folders, file names -------------------
        # experimentFolder = r'O:\HybridDune experiment\data RBR, OSSI\copy RBR Udrive series1'                  # path where the data is sitting # Rubens Laptop
        # subfolder_in_all = ['refP1 RBR4', 'S1P3 RBR5', 'S2P3 RBR1', 'S3P3 RBR6', 'S4P3 RBR2','S1P2 RBR3'] # subfolder where file is sitting within experimentFolder (on O drive Daan)
        # instrumentname_all = subfolder_in_all #['refP1 RBR4', 'S3P6 RBR6']
        # i = 4 #3 geeft error
        # instrument = instrumentname_all[i]  # select instrument
    ncInFile = os.path.join(input_folder, input_file)
    ncOutFile = os.path.join(output_folder,output_file)

    # Physical constants -------------------
    rho = 1027 #kg/m3, for seawater at 9C, avg temp at HKZ measurement station
    g = 9.8125  # value Zandmotor

    # Settings spectral analysis: segments (Welch method) -------------------
    p_blocks = 20          # number of segments within block, for the Welch method
    D_length = 1200        # Duration of block in seconds (20 minutes)
    D_length_s = '1200s'   # Duration of block in seconds (20 minutes). (same as above, but as string for xarray)


    # IMPORT DATA, RESAMPLE TO 2D, FILTER ===============================================================================================
    # NB: 3 datasets are created: ds0 (imported data), ds_2D (data reshaped to blocks of (e.g.) 20 minutes, and some temporary variables 
    # that won't be saved), and ds_out (final output dataset, with wave statistics and 1D timeseries of water level, depth)

    # load the data from netcdf ---------------------------------------------------
    ds0 = xr.open_dataset(ncInFile)    # dataset with relative pressure

    # frequency resolution in fourier space --------------------------------------------
    p_blocks = 20                                       # number of segments within block, for the Welch method
    D_length = 1200                                     # Duration of block in seconds (20 minutes)
    fresolution = p_blocks / D_length                   # Frequency resolution is 1/T_segment = n_segments / T_block
    nperseg = D_length * ds0.sf.values / p_blocks - 0.5 # dim should be len(ds.f); whelch has (nperseg/2 +1)

    # Make 2D data array with one column per burst -----------------------------------
    pt = ds0.p_rel.values                               # relative pressure, pAir subtracted
    nSamples = len(pt)
    dt = ds0.isel(t=1).t - ds0.isel(t=0).t

    burstDuration = pd.Timedelta(D_length_s)            # Burst duration (1200 seconds = 20 minutes)
    burstLength = int(burstDuration / dt)
    nBursts = int(np.floor(nSamples / burstLength))

    pt = pt[:nBursts * burstLength]
    t_full = ds0.t.values[:nBursts * burstLength]       # time vector for all samples, up to the last complete burst. skip incomplete burst at end
    t_block = t_full[::burstLength]                     # take every nth step, so t = t0 of every burst

    N = (ds0.t.values[:burstLength] - ds0.t.values[0]) / np.timedelta64(1, 's')  # New dimension: time in seconds since start of burst

    # Cast pressure into a 2D array 
    ds_2D = xr.Dataset(data_vars={},                    # Temporary 2D dataset, with cooridnates t (no. of blocks), N (obs within block)
                    coords={'t_full': t_full,                
                            't_block': t_block,
                            'N': N,
                            'f': np.arange(0, ds0.sf.values / 2, fresolution)})
    ds_2D['p'] = (('t_block', 'N'), pt.reshape((nBursts, burstLength)))      # relative pressure, pAir subtracted


    # make a new dataset that has an extra dimension to accomodate for the frequency axis ----------------------------------
    ds_out = xr.Dataset(data_vars={},
                    coords={'t_full': t_full,                  
                            't_block': t_block,
                            'f': np.arange(0, ds0.sf.values / 2, fresolution)})
    ds_out['t_full'].attrs = {'long name': 'time, of every observation'}
    ds_out['t_block'].attrs = {'long name': 'time of the start of every block (used for calculating wave statistics)'}
    ds_out['f'].attrs = {'units': 'Hz', 'long name': 'frequency'}


    # # put all variables in this new dataset
    ds0 = ds0.sel(t=slice(t_full[-1]))  # Make copy of ds0, with only complete blocks (i.e., 41 minutes data, 20min blocks, than take first 40 min)
    ds_out.attrs = ds0.attrs
    for key in ds0.data_vars:
        ds_out[key] = ds0[key]
    ds_out.attrs['comment_2'] = ds_out.attrs['comment_1'] # Made with XArray
    ds_out.attrs['comment_1'] = 'subscripts IG and WW in variable names refer to wave statistiscs of infragravity waves (f<0.05 Hz) or wind waves (f=0.05-1.5 Hz)'

    # Drop raw pressure, rename p_rel to p
    ds_out['p'] = ('t_full',ds_out.p_rel.data)          # renaming p_rel to p (overwriting existing p variable)
    ds_out['p'].attrs = ds0['p_rel'].attrs              # copying metadata
    ds_out = ds_out.drop_vars('p_rel')                  # and  dropping the old p_rel variable
    ds_out = ds_out.drop_dims('t')                      # plus dropping dimension t (now renamed to t_full)






    # Save to netcdf ====================================================================================================================
    #ds0.p.values = np.round(ds0.p.values)    # Round pressure to 1 Pa = 0.1 mm (smaller file size)
    #ds0.p_rel.values = np.round(ds0.p_rel.values)    
    encoding = {var: {"zlib": True, "complevel": 4} for var in list(ds_out.data_vars) + list(ds_out.coords)}  # Apply deflate compression to all variables and coordinates in netCDF

    # Save with compression
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    ds_out.to_netcdf(ncOutFile, encoding=encoding)

# ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    RUN WAVE PROCESSING    ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
wave_statistics_netcdf(input_folder, input_file, output_folder, output_file)
print(f"Processed data saved to {output_file}")
