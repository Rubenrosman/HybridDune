# read NETCDF file and plot the data
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# read NETCDF file
nc = Dataset('test.nc', 'r')
print(nc.variables.keys())

# read data
data = nc.variables['data'][:]
print(data.shape)

# plot data
plt.imshow(data)
plt.colorbar()
plt.show()

# close file
nc.close()

