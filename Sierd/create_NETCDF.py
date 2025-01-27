#create a NETCDF file with a variable and a dimension
import netCDF4
import numpy as np

# variables are time and cross shore distance
time = np.arange(0, 10, 1)
distance = np.arange(0, 100, 1)

# create a NETCDF file
nc = netCDF4.Dataset('test.nc', 'w')

# create a dimension
nc.createDimension('time', len(time))
nc.createDimension('distance', len(distance))

# create a variable
time_var = nc.createVariable('time', 'f4', ('time',))
distance_var = nc.createVariable('distance', 'f4', ('distance',))
data_var = nc.createVariable('data', 'f4', ('time', 'distance'))

#add x, y and z variables
x_var = nc.createVariable('x', 'f4', ('distance',))
y_var = nc.createVariable('y', 'f4', ('distance',))
z_var = nc.createVariable('z', 'f4', ('time','distance',))


# assign values to the variables
time_var[:] = time
distance_var[:] = distance
data_var[:] = np.random.rand(len(time), len(distance))
x_var[:] = np.random.rand(len(distance))
y_var[:] = np.random.rand(len(distance))
z_var[:] = np.random.rand(len(distance))



# close the file
nc.close()


