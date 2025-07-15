import os
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import xarray as xr
from scipy.ndimage import uniform_filter1d

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