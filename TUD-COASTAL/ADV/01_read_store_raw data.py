import os
import sys
sys.path.append(r'c:\checkouts\python\TUD-COASTAL\instrumentProcessing')
from vector import Vector
from datetime import datetime


# location of raw data
dataFolder = r'c:\checkouts\python\TUD-COASTAL\instrumentProcessing\example_data\ADV\raw_phzd'
# name of the instantiated vector class
name = 'vec1'
# start time over which to read data (must be larger than first recorded time)
tstart = '2020-11-30 17:00:00'
# stop time over which to read data (must be smaller than last recorded time)
tstop = '2020-12-02 17:00:00'

ncOutDir = r'c:\checkouts\python\TUD-COASTAL\instrumentProcessing\example_data\ADV\raw_netcdf'

# raw data to netcdf
vec = Vector(name, dataFolder, tstart=tstart, tstop=tstop)

# reads the raw data from tstart to tstop and casts all data in a pandas DataFrame that is stored under vec.dfpuv.
# in case there is no data between tstart and tstop the DataFrame is not instantiated
vec.read_raw_data()

# break up the data into burst blocks
vec.cast_to_blocks_in_xarray(blockWidth=600)

# compute burst averages (make sure to read vector.py what is happening exactly!)
vec.compute_block_averages()

# all data is collected in an xarray Dataset ds. We extract this from the class instantiation and
# we can easily write it to netCDF
ds = vec.ds

# add global attribute metadata
ds.attrs = {'Conventions': 'CF-1.6',
            'title': '{}'.format(vec.name),
            'instrument': '{}'.format('vec1'),
            'instrument serial number': '{}'.format(16725),
            'epsg': 28992,
             'x': 117196.6,
            'y': 559818.2,
            'time zone': 'UTC+2',
            'coordinate type': 'XYZ',
            'summary': 'December pilot field campaign',
            'contact person': 'Marlies van der Lugt',
            'emailadres': 'm.a.vanderlugt@tudelft.nl',
            'construction datetime': datetime.now().strftime("%d-%b-%Y (%H:%M:%S)"),
            'version': 'v1',
            'version comments': 'constructed with xarray'}

#specify compression for all the variables to reduce file size
comp = dict(zlib=True, complevel=5)
ds.encoding = {var: comp for var in ds.data_vars}

# save to netCDF
if not os.path.exists(ncOutDir):
    os.mkdir(ncOutDir)
ds.to_netcdf(ncOutDir + r'\{}_pilot.nc'.format(vec.name))
#








