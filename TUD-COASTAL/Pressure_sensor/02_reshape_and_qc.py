# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 08:27:20 2021

@author: marliesvanderl
"""
import os
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(r'c:\checkouts\python\TUD-COASTAL\python')
from KNMI_readers import read_knmi_uurgeg

experimentFolder = r'c:\checkouts\python\TUD-COASTAL\python\instrumentProcessing\example_data\SOLO'
rho = 1028
g = 9.8
zb = -0.5
zi =  -0.38
knmiFile = r'c:\Users\marliesvanderl\phd\vakbegeleiding\CIE5318\2022\SOLO\data\raw\uurgeg_235_2021-2030_dekooij.txt'

dataFile = os.path.join(experimentFolder,'raw_netcdf','L2C2SOLO.nc')
ds0 = xr.open_dataset(dataFile)
instr = ds0.instrument

# correct for the air pressure fluctuations and drift in the instrument
# first we load the data and add it to the dataset
dfp = read_knmi_uurgeg(knmiFile, stationNumber=235)
dt = ((ds0.t[1] - ds0.t[0]) / np.timedelta64(1, 's')).values # target frequency
pAir = dfp['P'].to_xarray().sel(t=slice(ds0.t.min(), ds0.t.max())).resample({'t': '{}S'.format(dt)}).interpolate('linear')
ds0['pAir'] = pAir.sel(t=slice(ds0.t.min(), ds0.t.max()))

# we correct for drift in air pressure, nothing else
ds0['dpAir'] = ds0['pAir'] - ds0['pAir'].isel(t=0)

# correct the pressure signal with dpAir and with drift in instrument pressure
ds0['pc'] = ds0['p'] -ds0['p'].min() - ds0['dpAir']
ds0['pc'].attrs = {'units': 'Pa + NAP', 'long_name': 'pressure', 'comments': 'drift in air pressure is corrected'}

# -----------------------------------------------------------------------------
# reshape to one row per burst in data array
pt = ds0.pc.values
nSamples = len(pt)
dt = ds0.isel(t=1).t - ds0.isel(t=0).t
sf = np.timedelta64(1, 's') / dt.values

burstDuration = pd.Timedelta('600S')
burstLength = int(burstDuration / dt)
nBursts = int(np.floor(nSamples / burstLength))

pt = pt[:nBursts * burstLength]
t = ds0.t[::burstLength]
t = t[:nBursts]
N = (ds0.t.values[:burstLength] - ds0.t.values[0]) / np.timedelta64(1, 's')
# pdb.set_trace()

# --------------------------------------------------------------------------
# cast into a 2D array
ds = xr.Dataset(data_vars={},
                coords={'t': t, 'N': N})
# copy all data over into this new structure
ds['p'] = (('t', 'N'), pt.reshape((nBursts, burstLength)))
ds['zi'] = zi
ds['zb'] = zb
ds['sf'] = sf


# remove all bursts where instrument fell dry
ds['p'] = ds.p.where(ds.p.std(dim='N') > 70)

# --------------------------------------------------------------------------
# pdb.set_trace()
ds['p'].attrs = {'units': 'Pa +NAP', 'long_name': 'pressure', 'comments': 'corrected for air pressure'}
ds['zi'].attrs = {'units': 'm+NAP', 'long_name': 'zi'}
ds['zb'].attrs = {'units': 'm+NAP', 'long_name': 'zb'}
ds['sf'].attrs = {'units': 'Hz', 'long_name': 'sampling frequency'}
ds.attrs = ds0.attrs
ds.attrs['summary'] = 'SEDMEX field campaign, pressure corrected for air pressure and cast in bursts of 10 minutes'
ds['name'] = instr
if not os.path.isdir(os.path.join(experimentFolder,'QC')):
    os.mkdir(os.path.join(experimentFolder,'QC'))
ncFilePath = os.path.join(experimentFolder, 'QC', instr + '.nc')
ds.to_netcdf(ncFilePath)

