# -*- coding: utf-8 -*-
"""
Created on Tue May 11 10:00:10 2021

@author: marliesvanderl
"""

import sys
sys.path.append(r'c:\checkouts\python\TUD-COASTAL\python\instrumentProcessing')
import os
import puv
import numpy as np
import xarray as xr
 
# input specification
instrFile = r'c:\checkouts\python\TUD-COASTAL\python\instrumentProcessing\example_data\ADV\qc\vec1.nc'
ncOutFile = r'c:\checkouts\python\TUD-COASTAL\python\instrumentProcessing\example_data\ADV\tailored\vec1_pilot_tailored.nc'
 
# frequency resolution in fourier space
fresolution = 0.03125
#number of directional bins 
ntheta = 64
  
# load the raw data from netcdf
ds0 = xr.open_dataset(instrFile).load()


# interpolate nans
N = len(ds0.N)
for var in ['u', 'v', 'p', 'eta']:
    # interpolate the bursts where there is less than 5% nans
    data = ds0[var].where(
        np.isnan(ds0[var]).sum(dim='N') < 0.05 * len(ds0.N)
    ).dropna(dim='t', how='all')
    if len(data.t) != 0:
        ds0[var] = data.interpolate_na(
            dim='N',
            method='cubic',
            max_gap=8)

    # and fill the gaps more than 8 in length with the burst average
    ds0[var] = ds0[var].fillna(ds0[var].mean(dim='N'))

ds0 = ds0.dropna(dim='t')

# make a new dataset that has an extra dimension to accomodate for the frequency axis
ds = xr.Dataset(data_vars={},
          coords = {'t': ds0.t.values,
                    'N': ds0.N.values,
                    'f': np.arange(0, ds0.sf.values/2, fresolution),
                    'theta': np.arange(start=-np.pi,stop=np.pi,step=2*np.pi/ntheta)})
ds['f'].attrs = {'units': 'Hz'}
ds.attrs = ds0.attrs

# put all variables in this new dataset
for key in ds0.data_vars:
    ds[key] = ds0[key]


# extract sampling frequency as explicit variable
sf = ds.f.values          

# compute water depth
ds['h'] = ds['zsmean']-ds['zb']   


# do several wave statistics computations, only based on pressure
ufunc = lambda x, h: puv.attenuation_corrected_wave_spectrum(
    'pressure',
    ds.sf.values, x, h,
    ds.zi.values,
    ds.zb.values,
    fresolution=fresolution)

fx, ds['vy'] = xr.apply_ufunc(ufunc,
                        ds['p'], ds['h'],
                        input_core_dims=[['N'], []],
                        output_core_dims=[['f'], ['f']],
                        vectorize=True) 
ds['vy'].attrs = {'units': 'm2/Hz', 'long_name': 'spectral density'}
               

ufunc = lambda vy: puv.get_peak_frequency(ds.f.values, vy)
ds['fp'] = xr.apply_ufunc(ufunc,
                        ds['vy'],
                        input_core_dims=[['f']],
                        output_core_dims=[[]], 
                        vectorize=True) 
 
ufunc = lambda vy, fp: puv.compute_wave_params(ds.f.values, vy, fmin=0.5*fp, fmax=5)
ds['Hm0'], ds['Tp'], ds['Tm01'], ds['Tm02'], ds['Tmm10'], ds['Tps'] = xr.apply_ufunc(ufunc,
                        ds['vy'], ds['fp'],
                        input_core_dims=[['f'], []],
                        output_core_dims=[[], [], [], [], [], []],
                        vectorize=True) 
ds['Hm0'].attrs = {'units': 'm', 'long_name': 'significant wave height','computation':'computed between fmin=0.5fp and fmax=5'}
ds['Tp'].attrs = {'units': 's', 'long_name': 'peak wave period','computation':'computed between fmin=0.5fp and fmax=5'}
ds['Tm01'].attrs = {'units': 's', 'long_name': 'mean wave period','computation':'computed between fmin=0.5fp and fmax=5'}
ds['Tm02'].attrs = {'units': 's', 'long_name': 'mean wave period','computation':'computed between fmin=0.5fp and fmax=5'}
ds['Tmm10'].attrs = {'units': 's', 'long_name': 'mean wave period','computation':'computed between fmin=0.5fp and fmax=5'}
ds['Tps'].attrs = {'units': 's', 'long_name': 'peak wave period','computation':'computed between fmin=0.5fp and fmax=5', 'comment':'smoothed estimate from the discrete spectrum'}


# compute current magnitudes and direction all computed in the time domain
ds['u_mean'] = ds.u.mean(axis=1)
ds['u_mean'].attrs = {'units': 'm/s', 'long_name': 'current x-component', 'computation': 'burst averaged'}

ds['v_mean'] = ds.v.mean(axis=1)
ds['v_mean'].attrs = {'units': 'm/s', 'long_name': 'current y-component', 'computation': 'burst averaged'}
               
ds['cur_dir'] = np.arctan2(ds['v_mean'], ds['u_mean'])*180/np.pi
ds['v_mean'].attrs = {'units': 'deg', 'long_name': 'current direction, cartesian convention'}

# directional wave spectra

ufunc = lambda p, u, v, h, fp: puv.wave_MEMpuv(p/1e4, u, v, h,
                    ds.zi.values,
                    ds.zb.values,
                    ds.sf.values,
                    fresolution=fresolution,
                    ntheta=ntheta,
                    fcorrmin=0.5*fp,
                    fcorrmax=5,
                    maxiter=20)
            
fx, vy, theta, ds['S'] = xr.apply_ufunc(ufunc,
                        ds['p'], ds['u'], ds['v'], ds['h'], ds['fp'],
                        input_core_dims=[['N'], ['N'], ['N'], [], []],
                        output_core_dims=[['f'], ['f'], ['theta'], ['f', 'theta']],
                        vectorize=True) 
ds['S'].attrs = {'units': 'm2/Hz/rad', 'long_name': 'directional variance density',
                 'computation': 'computed between fmin=0.5fp and fmax=5'}

# statistics from directional wave spectra
ufunc = lambda vy,S,fp: puv.compute_wave_params(ds.f.values, vy, fmin=0.5*fp, fmax=5, theta=ds.theta.values, S=S)
Hm0, Tp, Tm01, Tm02, Tmm10, Tps, ds['wavedirmean'],ds['dirspread'] = xr.apply_ufunc(ufunc,
                        ds['vy'], ds['S'], ds['fp'],
                        input_core_dims=[['f'], ['f', 'theta'], []],
                        output_core_dims=[[], [], [], [], [], [], [], []],
                        vectorize=True) 
ds['wavedirmean'].attrs = {'units': 'deg', 'long_name': 'mean wave direction', 'computation': 'computed between fmin=0.5fp and fmax=5'}
ds['dirspread'].attrs = {'units': 'deg', 'long_name': 'directional spreading', 'computation': 'computed between fmin=0.5fp and fmax=5'}
           
        
# write to file
# we strip all information on burst scale from the dataset to reduce size (and this info is already present in the raw_netcdf version of the data)
dsTailored = ds.drop_dims('N')
dsTailored.to_netcdf(ncOutFile)

 
