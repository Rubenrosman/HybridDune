# ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    Import modules    ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    
# This script processes wave data from a netCDF file, calculates wave statistics and saves them to a new netCDF fifle.
# First the fuction for wave processing is defined, then the function is called with the input and output file names.

import os
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import xarray as xr
import matplotlib.pyplot as plt
import sys
from scipy.signal import welch
#sys.path.append(r'C:\Users\dpoppema\Documents\GitHub\HybridDune\Ruben\Pressure_sensors\S1\RBR_05') # to find the puv.py file
import puv 


# ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    DEFINE FUNCTION FOR WAVE PROCESSING   ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
def wave_statistics_netcdf(input_folder, input_file, output_folder, output_file):
    """
    Function to process wave data from a netCDF file, resample it to 2D blocks, and calculate wave statistics.
    The function reads the input data, reshapes it, calculates wave statistics using puv.py and saves the processed data to a new netCDF file.
    """
    # ### INPUT/OUTPUT FILES, SETTINGS ################################################################################################################################
    # Folders, file names -------------------
    ncInFile = os.path.join(input_folder, input_file)
    ncOutFile = os.path.join(output_folder,output_file)

    # Physical constants -------------------
    rho = 1027 #kg/m3, for seawater at 9C, avg temp at HKZ measurement station
    g = 9.8125  # value Zandmotor

    # Settings spectral analysis: segments (Welch method) -------------------
    p_blocks = 20          # number of segments within block, for the Welch method
    D_length = 1200        # Duration of block in seconds (20 minutes)
    D_length_s = '1200s'   # Duration of block in seconds (20 minutes). (same as above, but as string for xarray)


    # ### IMPORT DATA, RESAMPLE TO 2D  ################################################################################################################################
    # NB: 3 datasets are created: ds0 (imported data), ds_2D (data reshaped to blocks of (e.g.) 20 minutes, and some temporary variables 
    # that won't be saved), and ds_out (final output dataset, with wave statistics and 1D timeseries of water level, depth)

    # load the data from netcdf ---------------------------------------------------
    ds0 = xr.open_dataset(ncInFile)    # dataset with relative pressure

    # frequency resolution in fourier space --------------------------------------------
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
    ds_out.attrs['summary'] = 'Hybrid-Dune campaign: pressure (corrected for air pressure) and wave statistics)'

    # Drop raw pressure, rename p_rel to p
    ds_out['p'] = ('t_full',ds_out.p_rel.data)          # renaming p_rel to p (overwriting existing p variable)
    ds_out['p'].attrs = ds0['p_rel'].attrs              # copying metadata
    ds_out = ds_out.drop_vars('p_rel')                  # and  dropping the old p_rel variable
    ds_out = ds_out.drop_dims('t')                      # plus dropping dimension t (now renamed to t_full)


    # ### APPLY BANDPASS FILTER, DEAL WITH TIME-VARYING BED/INSTRUMENT LEVELS AND COMPUTE BLOCK-AVERAGED WATER DEPTH ################################################################
    # pressure: bandpass filtering  --------------------------------------------------------
    ufunc_bandpass    = lambda x: puv.band_pass_filter2(ds0.sf.values, x, fmin=0.004, fmax=1.5)
    ufunc_bandpass_IG = lambda x: puv.band_pass_filter2(ds0.sf.values, x, fmin=0.004, fmax=0.05)
    ufunc_bandpass_WW = lambda x: puv.band_pass_filter2(ds0.sf.values, x, fmin=0.05, fmax=1.5)

    ds_2D['p'] = xr.apply_ufunc(ufunc_bandpass, 
                            ds_2D['p'],
                            input_core_dims=[['N']],
                            output_core_dims=[['N']],
                            vectorize=True)
    ds_2D['p_IG'] = xr.apply_ufunc(ufunc_bandpass_IG, 
                            ds_2D['p'],
                            input_core_dims=[['N']],
                            output_core_dims=[['N']],
                            vectorize=True)
    ds_2D['p_WW'] = xr.apply_ufunc(ufunc_bandpass_WW, 
                            ds_2D['p'],
                            input_core_dims=[['N']],
                            output_core_dims=[['N']],
                            vectorize=True)

    # # Time-varying bed and instrument level --------------------------------------------------------
    # Define z_bed_obs and t_bed_obs for RBR6 
    t_zb_obs = pd.to_datetime(ds0.t_zb.values)
    t_block2 = pd.to_datetime(t_block)  # Convert to pandas datetime for consistency
    zb_obs = ds0.zb.values

    # Interpolate z_bed_obs to z_bed_block for time t_block
    zb_block = np.interp(t_block2, t_zb_obs, zb_obs)            # Interpolate bed level to block time vector
    zb_block2 = np.reshape(zb_block,(len(zb_block),1))          # reshape from row to column vector for compatibility with p

    # repeat for z_i. Only for file no. 3 (RBR4), the rest had constant bed level, so simply repeat bedlevel there
    if np.size(ds0.zi.values) > 1: # If the instrument was moved, so time-varying zi
        t_zi_obs = pd.to_datetime(ds0.t_zi.values)
        zi_obs = ds0.zi.values

        # 'Interpolate' zi to t_block: use previous known value (so step-wise, not linear interpolation)
        idx = np.searchsorted(t_zi_obs, t_block2, side='right') - 1
        idx = np.clip(idx, 0, len(zi_obs) - 1)                  # Ensure indices are within bounds
        zi_block = zi_obs[idx]
        zi_block2 = np.reshape(zi_block, (len(zi_block), 1))    # reshape from row to column vector for compatibility with p
    else: # if zi did not vary during test, simply repeat zi
        zi_block = np.repeat(ds0.zi.values, len(zb_block))
        zi_block2 = np.reshape(zi_block, (len(zi_block), 1))    # reshape from row to column vector for compatibility with p

    # Filter pressure for negative pressures (filter 1), then calculate water depth ----------------------------------------------------
    block_mask = ds_2D['p'] < 0
    ds_2D['p'] = ds_2D['p'].where(~block_mask, 0)               # set negative pressures to zero

    # Compute water depth ----------------------------------------------------
    ds_2D['h_mean'] = ( (ds_2D['p'] / rho / g) + zi_block2 - zb_block2 ).mean(dim='N') # Mean water depth per burst: h=p/rho/g + (z_i-z_b)
    ds_out['h_mean'] = ds_2D['h_mean']
    ds_out['h_mean'].attrs = {'long_name': 'mean water depth', 'units': '[m]'} # avg per burst


    # ### COMPUTE WAVE STATISTICS. 3 TIMES: FOR FULL BANDPASS, INFRAGRAFITY AND WIND WAVES ###############################################################################
    # Attenuate signal: from pressure to surface elevation -------------------------------------
    ufunc_attenuate = lambda x, h, zi, zb: puv.attenuate_signal(
        'pressure', 
        ds0.sf.values, x, h, 
        zi,
        zb, 
        rho = rho,
        g = g,
        removeNoise=False,
        detrend=True)

    fx, h = xr.apply_ufunc(ufunc_attenuate,                                        # attenuate signal (calculate water surface), using p(t), and h_mean and zb_mean (each 1 value per block)
                                ds_2D['p'], ds_2D['h_mean'], zi_block, zb_block,
                                input_core_dims=[['N'], [],[],[]],
                                output_core_dims=[['f'], ['N']],
                                vectorize=True)
    fx, h_IG = xr.apply_ufunc(ufunc_attenuate,                                     # repeat for infragrafity waves
                                ds_2D['p_IG'], ds_2D['h_mean'], zi_block, zb_block,
                                input_core_dims=[['N'], [],[],[]],
                                output_core_dims=[['f'], ['N']],
                                vectorize=True)
    fx, h_WW = xr.apply_ufunc(ufunc_attenuate,                                     # repeat for wind waves
                                ds_2D['p_WW'], ds_2D['h_mean'], zi_block, zb_block,
                                input_core_dims=[['N'], [],[],[]],
                                output_core_dims=[['f'], ['N']],
                                vectorize=True)

    # Calculate water surface elevation from depth -------------------------------------
    ds_2D['zs']    = zb_block2 + h                                                 # water level (surface elevation) = bed level + depth
    ds_2D['zs_IG'] = zb_block2 + h_IG
    ds_2D['zs_WW'] = zb_block2 + h_WW

    ds_out['zs']    = ('t_full',np.ravel(ds_2D['zs'])*1000)                        # convert to mm, for smaller file size (rounding to mm, done later, gives integers)
    ds_out['zs_IG'] = ('t_full',np.ravel(ds_2D['zs_IG'])*1000)                        
    ds_out['zs_WW'] = ('t_full',np.ravel(ds_2D['zs_WW'])*1000)                        
    ds_out['h']     = ('t_full',np.ravel(h)*1000)                                  # water depth, converted to mm
    ds_out['h_IG']  = ('t_full',np.ravel(h_IG)*1000)
    ds_out['h_WW']  = ('t_full',np.ravel(h_WW)*1000)

    ds_out['zs'].attrs    = {'units': 'mm +NAP', 'long_name': 'surface elevation'}   # DAAN: CHECK: MORE METADATA NEEDED?
    ds_out['zs_IG'].attrs = {'units': 'mm +NAP', 'long_name': 'surface elevation infragrafity waves'}   
    ds_out['zs_WW'].attrs = {'units': 'mm +NAP', 'long_name': 'surface elevation wind waves'}   
    ds_out['h'].attrs     = {'units': 'mm', 'long_name': 'water depth'}              
    ds_out['h_IG'].attrs  = {'units': 'mm', 'long_name': 'water depth infragrafity waves'}
    ds_out['h_WW'].attrs  = {'units': 'mm', 'long_name': 'water depth wind waves'}

    # Filter 2: make blocks NaN if more than 25% of the observations had a depth_above_instrument < 0.05 m
    block_mask = ((ds_2D['zs'] - zi_block2) < 0.05).mean(dim='N') >= 0.25 # to be used for for Welch, calculation wave properties:
    z_filtered = ds_2D['zs']
    z_filtered_IG = ds_2D['zs_IG']
    z_filtered_WW = ds_2D['zs_WW']
    z_filtered    = z_filtered.where(~block_mask, np.nan)
    z_filtered_IG = z_filtered_IG.where(~block_mask, np.nan)
    z_filtered_WW = z_filtered_WW.where(~block_mask, np.nan)

    # Determine wave spectrum -------------------------------------
    ufunc_welch = lambda p: welch(p, fs=ds0.sf.values, nperseg=nperseg, detrend='constant', window='hann') # Detrend: is per segment. 1min, about constant, so false used

    ds_2D['frequencies'], ds_2D['psd'] = xr.apply_ufunc(ufunc_welch,
                                                    z_filtered,
                                                    input_core_dims=[['N']],
                                                    output_core_dims=[['f'], ['f']],
                                                    vectorize=True)
    ds_2D['frequencies'], ds_2D['psd_IG'] = xr.apply_ufunc(ufunc_welch,
                                                    z_filtered_IG,
                                                    input_core_dims=[['N']],
                                                    output_core_dims=[['f'], ['f']],
                                                    vectorize=True)
    ds_2D['frequencies'], ds_2D['psd_WW'] = xr.apply_ufunc(ufunc_welch,
                                                    z_filtered_WW,
                                                    input_core_dims=[['N']],
                                                    output_core_dims=[['f'], ['f']],
                                                    vectorize=True)
    ds_out['psd'] = ds_2D['psd']
    ds_out['psd_IG'] = ds_2D['psd_IG']
    ds_out['psd_WW'] = ds_2D['psd_WW']
    ds_out['psd'].attrs = {'units': 'm^2/Hz', 'long_name': 'Power spectral density (Welch) (f=0.004-1.5 Hz)'}
    ds_out['psd_IG'].attrs = {'units': 'm^2/Hz', 'long_name': 'Power spectral density (Welch) of long waves (<0.05 Hz)'}
    ds_out['psd_WW'].attrs = {'units': 'm^2/Hz', 'long_name': 'Power spectral density (Welch) of wind waves (0.05-1.5 Hz)'}

    # Determine wave parameters ---------------------------
    ufunc_wave_params = lambda psd: puv.compute_wave_params(ds_2D.f.values, psd, fmin=0.004 , fmax=5)   

    ds_out['Hm0'], ds_out['Tp'], ds_out['Tm01'], ds_out['Tm02'], ds_out['Tmm10'], ds_out['Tps'] = xr.apply_ufunc(ufunc_wave_params,
                                                                            ds_2D['psd'],
                                                                            input_core_dims=[['f']],
                                                                            output_core_dims=[[], [], [], [], [], []],
                                                                            vectorize=True)
    ds_out['Hm0_IG'], ds_out['Tp_IG'], ds_out['Tm01_IG'], ds_out['Tm02_IG'], ds_out['Tmm10_IG'], ds_out['Tps_IG'] = xr.apply_ufunc(ufunc_wave_params,
                                                                            ds_2D['psd_IG'],
                                                                            input_core_dims=[['f']],
                                                                            output_core_dims=[[], [], [], [], [], []],
                                                                            vectorize=True)
    ds_out['Hm0_WW'], ds_out['Tp_WW'], ds_out['Tm01_WW'], ds_out['Tm02_WW'], ds_out['Tmm10_WW'], ds_out['Tps_WW'] = xr.apply_ufunc(ufunc_wave_params,
                                                                            ds_2D['psd_WW'],
                                                                            input_core_dims=[['f']],
                                                                            output_core_dims=[[], [], [], [], [], []],
                                                                            vectorize=True)
    ds_out['Hm0'].attrs = ds_out['Hm0_IG'].attrs = ds_out['Hm0_WW'].attrs = {'units': 'm', 'long_name': 'significant wave height: Hm0=4sqrt(m0), with m0 zeroth-order spectral moment'}
    ds_out['Tp'].attrs = ds_out['Tp_IG'].attrs = ds_out['Tp_WW'].attrs = {'units': 's', 'long_name': 'peak wave period'}
    ds_out['Tm01'].attrs = ds_out['Tm01_IG'].attrs = ds_out['Tm01_WW'].attrs = {'units': 's', 'long_name': 'mean wave period: Tm01 = m0/m1'}
    ds_out['Tm02'].attrs = ds_out['Tm02_IG'].attrs = ds_out['Tm02_WW'].attrs = {'units': 's', 'long_name': 'mean zero-crossing period: Tm02 = sqrt(m0/m2)'}
    ds_out['Tmm10'].attrs = ds_out['Tmm10_IG'].attrs = ds_out['Tmm10_WW'].attrs = {'units': 's', 'long_name': 'T-1,0: mean absolute wave period T-1,0 = m_(-1)/m_0'}
    ds_out['Tps'].attrs = ds_out['Tps_IG'].attrs = ds_out['Tps_WW'].attrs = {'units': 's', 'long_name': 'Smoothed peak wave period'}


    # ### SKEWNESS AND ASYMMETRY ################################################################################################################################
    ## skewness of waves ##
    ufunc = lambda p: puv.compute_SkAs(ds0.sf.values,p,fpfac =None, fbounds = None)     # Daan: check calculation, metadata

    ds_out['Sk'], ds_out['As'], ds_out['sigma'] =  xr.apply_ufunc(ufunc,
                                                    ds_2D['p'], 
                                                    input_core_dims=[['N']],
                                                    output_core_dims=[[], [], []],
                                                    vectorize=True)
    ds_out['Sk'].attrs = {'units': 'm3', 'long name': 'wave skewness'}
    ds_out['As'].attrs = {'units': 'm3', 'long name': 'wave asymmetry'}
    ds_out['sigma'].attrs = {'units': 'm', 'long name': 'standard deviation'}


    # ### Save to netcdf ################################################################################################################################  
    # Compression settings
    ds_out.zs.values    = np.round(ds_out.zs.values)          # Save water level (in mm) with 1 mm accuracy. Decreases file size: less info, and datta can be saved as integers
    ds_out.zs_IG.values = np.round(ds_out.zs_IG.values)    
    ds_out.zs_WW.values = np.round(ds_out.zs_WW.values)    
    ds_out.h.values     = np.round(ds_out.h.values)           # Idem for water depth
    ds_out.h_IG.values  = np.round(ds_out.h_IG.values)    
    ds_out.h_WW.values  = np.round(ds_out.h_WW.values)

    encoding = {var: {"zlib": True, "complevel": 4} for var in list(ds_out.data_vars) + list(ds_out.coords)}  # Apply deflate compression to all variables and coordinates in netCDF

    # Save with compression
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    ds_out.to_netcdf(ncOutFile, encoding=encoding)


# ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    DEFINE INPUT/OUTPUT FILES AND RUN FUNCTION  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    
# RBR files, period 1 -------------------------------------------------------------------------------------------------------------------------------
input_folder = r'O:\HybridDune experiment\data RBR, OSSI\copy RBR Udrive series1\QC'
output_folder = r'O:\HybridDune experiment\data RBR, OSSI\copy RBR Udrive series1\processed'
names_all = ['S1P3 RBR5', 'S2P3 RBR1', 'S3P3 RBR6', 'S4P3 RBR2','S1P2 RBR3'] # part of file name that refers to sensor

for i in range(0,1):
    input_file = 'Pressure sensor ' + names_all[i] + ' p_rel - period 1.nc'
    output_file = 'Pressure sensor ' + names_all[i] + ' processed data - period 1 new.nc'

    # Call the function to process the wave data
    wave_statistics_netcdf(input_folder, input_file, output_folder, output_file)
    print(f"Processed data saved to {output_file}")


# # RBR files, period 2 -------------------------------------------------------------
# input_folder = r'O:\HybridDune experiment\data RBR, OSSI\copy RBR Udrive series2\QC'
# output_folder = r'O:\HybridDune experiment\data RBR, OSSI\copy RBR Udrive series2\processed'
# names_all = ['S3P3 RBR6', 'S1P2 RBR2', 'S3P2 RBR4', 'S4P2 RBR1','S3P1 RBR5'] # part of file name that refers to sensor

# for i in range(0,5):
#     input_file = names_all[i] + ' p_rel - period 2.nc'
#     output_file = 'Pressure sensor ' + names_all[i] + ' processed data - period 2.nc'

#     # Call the function to process the wave data
#     # wave_statistics_netcdf(input_folder, input_file, output_folder, output_file)
#     # print(f"Processed data saved to {output_file}")


# # OSSI files, period 1 -------------------------------------------------------------
# input_folder = r'O:\HybridDune experiment\data RBR, OSSI\Ossi data\QC'
# output_folder = r'O:\HybridDune experiment\data RBR, OSSI\Ossi data\processed'
# names_all = ['S2P2 Ossi1', 'S3P4 Ossi2', 'S4P2 Ossi3',  'S1P1 Ossi4', 'S2P1 Ossi5', 'S3P1 Ossi6', 'S4P1 Ossi7'] # part of file name that refers to sensor

# for i in range(0,7):
#     input_file = names_all[i] + ' p_rel - period 1.nc'
#     output_file = 'Pressure sensor ' + names_all[i] + ' processed data - period 1.nc'

#     # Call the function to process the wave data
#     wave_statistics_netcdf(input_folder, input_file, output_folder, output_file)
#     print(f"Processed data saved to {output_file}")
