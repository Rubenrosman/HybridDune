# -*- coding: utf-8 -*-
"""
Created on Tue May 11 10:00:10 2021

@author: marliesvanderl
"""

import sys
import os
sys.path.append(r'c:\checkouts\python\TUD-COASTAL\python')
import puv
import numpy as np
import xarray as xr

# %% input specification
experimentFolder = r'c:\checkouts\python\TUD-COASTAL\python\instrumentProcessing\example_data\SOLO'
instrFile = os.path.join(experimentFolder,'QC\L2C2SOLO.nc')
ncOutFile = os.path.join(experimentFolder,'tailored\L2C2SOLO.nc')

# frequency resolution in fourier space
fresolution = 0.03125
rho = 1025
g = 9.8

# %% load the raw data from netcdf
ds0 = xr.open_dataset(instrFile)

#let's remove the bursts where there are only nans
ds0 = ds0.dropna(dim='t')

# make a new dataset that has an extra dimension to accomodate for the frequency axis
ds = xr.Dataset(data_vars={},
                coords={'t': ds0.t.values,
                        'N': ds0.N.values,
                        'f': np.arange(0, ds0.sf.values / 2, fresolution)})
ds['f'].attrs = {'units': 'Hz'}
ds.attrs = ds0.attrs

# put all variables in this new dataset
for key in ds0.data_vars:
    ds[key] = ds0[key]

# extract sampling frequency as explicit variable
sf = ds.f.values

# compute water depth
ds['h'] = (ds['p']/rho/g + ds['zi']).mean(dim='N')
ds['h'].attrs = {'long_name': 'mean water level', 'units': 'm+NAP'}

# %% do several wave statistics computations, only based on pressure
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

ufunc = lambda vy, fp: puv.compute_wave_params(ds.f.values, vy, fmin=0.5 * fp, fmax=5)
ds['Hm0'], ds['Tp'], ds['Tm01'], ds['Tm02'], ds['Tmm10'], ds['Tps'] = xr.apply_ufunc(ufunc,
                                                                          ds['vy'], ds['fp'],
                                                                          input_core_dims=[['f'], []],
                                                                          output_core_dims=[[], [], [], [], [], []],
                                                                          vectorize=True)
ds['Hm0'].attrs = {'units': 'm', 'long_name': 'significant wave height',
                   'computation': 'computed between fmin=0.5fp and fmax=5'}
ds['Tp'].attrs = {'units': 's', 'long_name': 'peak wave period',
                  'computation': 'computed between fmin=0.5fp and fmax=5'}
ds['Tm01'].attrs = {'units': 's', 'long_name': 'mean wave period',
                    'computation': 'computed between fmin=0.5fp and fmax=5'}
ds['Tm02'].attrs = {'units': 's', 'long_name': 'mean wave period',
                    'computation': 'computed between fmin=0.5fp and fmax=5'}
ds['Tmm10'].attrs = {'units': 's', 'long_name': 'mean wave period',
                     'computation': 'computed between fmin=0.5fp and fmax=5'}
ds['Tps'].attrs = {'units': 's', 'long_name': 'peak wave period','computation':'computed between fmin=0.5fp and fmax=5', 'comment':'smoothed estimate from the discrete spectrum'}


# %% write to file
# we strip all information on burst scale from the dataset to reduce size (and this info is already present in the raw_netcdf version of the data)
dsTailored = ds.drop_dims('N')
if not os.path.isdir(os.path.join(experimentFolder,'tailored')):
    os.mkdir(os.path.join(experimentFolder,'tailored'))
ncFilePath = os.path.join(experimentFolder, 'tailored', ds0.instrument + '.nc')
dsTailored.to_netcdf(ncFilePath)

