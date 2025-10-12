import sys
sys.path.append(r'c:\checkouts\python\TUD-COASTAL\python\instrumentProcessing')
import os
import numpy as np
import xarray as xr
import puv
from KNMI_readers import read_knmi_uurgeg

# general settings

# location of knmi file (to correct for air pressure drift during the experiment)
knmiFile = r'c:\checkouts\python\TUD-COASTAL\python\instrumentProcessing\example_data\KNMI_20201208_hourly.txt'
# number of the knmi station (to make sure the correction is done with the correct KNMI station)
stationNumber = 235
# height of instrument above bed
hi = 0.57 #m
# height of instruments pressure sensor above bed
hip = 0.27 #m
# bed level
zb = -1.13#m NAP
# angle of x-pod of the vector head with respect to north (clockwise positive)
thet = 45
# density of water
rho = 1025 # kg/m3
# gravitational acceleration
g = 9.81 # m/s2
# parameters for the quality control:
QC = {
     'uLim':2.1, #maximum acceptable recorded u-velocity
     'vLim':2.1, #maximum acceptable recorded v-velocity
     'wLim':0.6, #maximum acceptable recorded w-velocity
     'corTreshold':70, #minimum correlation
     'maxFracNans': 0.02, #maximum fraction of rejected pings in the sample to proceed with processing based on interpolation
     'maxGap' : 4 #maximum amount of sequential rejected pings in the sample to proceed with processing based on interpolation
      }

ds = xr.open_dataset(r'c:\checkouts\python\TUD-COASTAL\python\instrumentProcessing\example_data\ADV\raw_netcdf\vec1_pilot.nc')

# % add some data to the dataset
ds['zb'] = zb

ds['zi'] = ds['zb'] + hi
ds['zi'].attrs = {'units': 'm+NAP', 'long_name': 'position probe'}

ds['zip'] = ds['zb'] + hip
ds['zip'].attrs = {'units': 'm+NAP', 'long_name': 'position pressure sensor'}

ds['rho'] = rho
ds['rho'].attrs = {'units': 'kg/m3', 'long_name': 'water density'}

ds['g'] = g
ds['g'].attrs = {'units': 'm', 'long_name': 'gravitational acceleration'}

# if correlation is outside confidence range
mc1 = ds.cor1 > QC['corTreshold']
mc2 = ds.cor2 > QC['corTreshold']
mc3 = ds.cor3 > QC['corTreshold']

# if observation is outside of velocity range
mu1 = np.abs(ds.u) < QC['uLim']
mu2 = np.abs(ds.v) < QC['uLim']
mu3 = np.abs(ds.w) < QC['uLim']

# if du larger than 4*std(u) then we consider it outlier and hence remove:
md1 = np.abs(ds.u.diff('N')) < 3 * ds.u.std(dim='N')
md1 = md1.combine_first(mu1)
md2 = np.abs(ds.v.diff('N')) < 3 * ds.v.std(dim='N')
md2 = md1.combine_first(mu2)
md3 = np.abs(ds.w.diff('N')) < 3 * ds.w.std(dim='N')
md3 = md1.combine_first(mu3)

ds['mc'] = np.logical_and(np.logical_and(mc1, mc2), mc3)
ds['mu'] = np.logical_and(np.logical_and(mu1, mu2), mu3)
ds['md'] = np.logical_and(np.logical_and(md1, md2), md3)
ds['mc'].attrs = {'units': '-', 'long_name': 'mask correlation'}
ds['mu'].attrs = {'units': '-', 'long_name': 'mask vel limit'}
ds['md'].attrs = {'units': '-', 'long_name': 'mask deviation'}

mp = np.abs(ds.p.diff('N')) < 4 * ds.p.std(dim='N')
mp = xr.concat([mp.isel(N=0), mp], dim="N")

ds.coords['maskp'] = (('t', 'N'), mp.values)
ds.coords['maskv'] = (('t', 'N'), np.logical_and(np.logical_and(ds.mc.values, ds.mu.values), ds.md.values))

# correct for the air pressure fluctuations and drift in the instrument
# first we load the data and add it to the dataset
dfp = read_knmi_uurgeg(
    knmiFile,
    stationNumber)
dt = ((ds.t[1] - ds.t[0]) / np.timedelta64(1, 's')).values
pAir = dfp['P'].to_xarray().resample({'t': '{}S'.format(dt)}).interpolate('linear')
ds['pAir'] = pAir.sel(t=slice(ds.t.min(), ds.t.max()))

# we correct for drift in air pressure, nothing else
ds['dpAir'] = ds['pAir'] - ds['pAir'].isel(t=0)

# correct the pressure signal with dpAir and with drift in instrument pressure
ds['pc'] = ds['p'] - ds['dpAir']
ds['pc'].attrs = {'units': 'Pa + NAP', 'long_name': 'pressure', 'comments': 'drift in air pressure is corrected'}

ds['eta'] = ds['pc'] / rho / g + ds.zip
ds['eta'].attrs = {'units': 'm+NAP', 'long_name': 'hydrostatic water level'}

ds['zsmean'] = ds.eta.mean(dim='N')
ds['zsmean'].attrs = {'units': 'm + NAP', 'long_name': 'water level',
                      'comments': 'burst averaged'}

ds['h'] = ds.zsmean - zb
ds['h'].attrs = {'units': 'm', 'long_name': 'water column height'}

# #% rotate to ENU coordinates (this is only necessary if measurements were performed in XYZ
ufunc = lambda u,v: puv.rotate_velocities(u,v,thet-90)
ds['u'],ds['v'] = xr.apply_ufunc(ufunc,
                    ds['u'], ds['v'],
                    input_core_dims=[['N'], ['N']],
                    output_core_dims=[['N'],['N']],
                    vectorize=True)
ds['u'].attrs = {'units':'m/s','long_name':'velocity E'}
ds['v'].attrs = {'units':'m/s','long_name':'velocity N'}
ds['w'].attrs = {'units':'m/s','long_name':'velocity U'}

# remove pressure observations where the estimated water level is
# lower than the sensor height with margin of error of 10 cm
ds.coords['maskd'] = (('t', 'N'), zb+hi < (ds['eta'].values - 0.1))
ds[['u','v','w','p','pc','eta']] = ds[['u','v','w','p','pc','eta']].where(ds.maskp == True)
ds[['u','v','w','p','pc','eta']] = ds[['u','v','w','p','pc','eta']].where(ds.maskd == True)
ds[['u','v','w','p','pc','eta']] = ds[['u','v','w','p','pc','eta']].where(ds.maskv == True)

# ammending the meta data to add extra info
ds.attrs['version'] = 'v2'
ds.attrs['coordinate type'] = 'ENU'
ds.attrs['comment'] = 'Quality checked data: pressure reference level corrected for airpressure drift,' + \
                 r'correlation and amplitude checks done and spikes were removed. ' + \
                 r'Velocities rotated to ENU coordinates based on heading and configuration in the field.'

# save to netCDF wwhere we don't include the sen data any more because we have only used it for the quality check
ds = ds.drop(['a1', 'a2', 'a3',
              'cor1', 'cor2', 'cor3',
              'snr1', 'snr2', 'snr3',
              'heading', 'pitch', 'roll',
              'voltage', 'pc'])

# specify compression for all the variables to reduce file size
comp = dict(zlib=True, complevel=5)
ds.encoding = {var: comp for var in ds.data_vars}
ncOutDir = r'c:\checkouts\python\TUD-COASTAL\python\instrumentProcessing\example_data\ADV\qc'
if not os.path.exists(ncOutDir):
    os.mkdir(ncOutDir)
ds.to_netcdf(os.path.join(ncOutDir, 'vec1.nc'), encoding=ds.encoding)

