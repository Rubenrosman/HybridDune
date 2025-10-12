import os
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import xarray as xr
from scipy.ndimage import uniform_filter1d
import sys
sys.path.append(r'C:\Users\dpoppema\Documents\GitHub\HybridDune\Daan\Python\RBR')
import wave_netcdf_func


fruits = ["apple", "banana", "cherry"]
for x in fruits:
  print(x)

a = [1, 1, 1, 2, 2, 3, 4]

#def movmean(x,window):
#    avg = np.convolve(x, np.ones(window)/window, mode='same')
#    return avg
#print(movmean(a,3))

y = uniform_filter1d(a, size=3)
print(y)

input_folder = r'O:\HybridDune experiment\data RBR, OSSI\copy RBR Udrive series1\QC'
input_file = 'S4P3 RBR2 p_rel - period 1.nc'
output_folder = r'O:\HybridDune experiment\data RBR, OSSI\copy RBR Udrive series1\processed'
output_file = 'Pressure sensor S4P3 RBR2 processed data - period 1.nc'

wave_netcdf_func.wave_statistics_netcdf(input_folder, input_file, output_folder, output_file)
print(f"Processed data saved to {output_file}")
