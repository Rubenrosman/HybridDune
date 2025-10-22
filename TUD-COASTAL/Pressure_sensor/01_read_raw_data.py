from datetime import datetime
import os
import numpy as np
import pandas as pd

def solo_data_reader(dataFile, sf):
    '''
    Function to read solo datafile.
    Returns a dataframe with a time column and pressure column in Pascal
    '''
    p = []
    datt = []
    with open(dataFile) as myfile:
        for index, line in enumerate(myfile):
            if index >= 1:
                lin = line.split(',')
                datt.append(lin[0])
                p.append(float(lin[1]))
    p = np.array(p) * 1e4  # dBar to Pa

    t = pd.date_range(datt[0], periods=len(datt), freq='{}S'.format(1 / sf))

    dfp = pd.DataFrame(data={'p': p}, index=t)

    dfp.index.name = 't'
    return dfp

#######################################################################################################################
# input parameters
sf = 8  # sampling frequency
experimentFolder = r'c:\checkouts\python\TUD-COASTAL\python\instrumentProcessing\example_data\SOLO'  # path where the data is sitting
instrumentName = 'L2C2SOLO'  # designated name of the instrument
dataFile =  os.path.join(experimentFolder, r'202440_20210919_1930_data.txt') # name of the datafile
serial_number = '18.09.00.08'  # unique serial number of the instrument
xRD = 117221.8  # x position of placement in field
yRD =  559793.1 # y position of placement in field

#do the reading from file and cast in xarray dataset
dfp = solo_data_reader(dataFile,sf)
ds = dfp.to_xarray()
ds.p.attrs = {'long_name': 'pressure', 'units': 'Pa'}

# add global attribute metadata
ds.attrs = {
    'Conventions': 'CF-1.6',
    'name': '{}'.format(instrumentName),
    'instrument': '{}'.format(instrumentName),
    'instrument type': 'Ruskin RBR Solo',
    'instrument serial number': '{}'.format(serial_number),
    'epsg': 28992,
    'x': xRD,
    'y': yRD,
    'time zone': 'UTC+2',
    'coordinate type': 'XYZ',
    'summary': 'SEDMEX field campaign',
    'contact person': 'Marlies van der Lugt',
    'emailadres': 'm.a.vanderlugt@tudelft.nl',
    'construction datetime': datetime.now().strftime("%d-%b-%Y (%H:%M:%S)"),
    'version': 'v1',
    'version comments': 'constructed with xarray'}

# save to netcdf
ncOutDir = os.path.join(experimentFolder,'raw_netcdf')
if not os.path.isdir(ncOutDir):
    os.mkdir(ncOutDir )
ds.to_netcdf(os.path.join(ncOutDir, instrumentName + '.nc'))