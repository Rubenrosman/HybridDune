# Script for importing netcdf, finding optimal shuffle settings for encoding, and then saving the dataset with these settings
import xarray as xr
import os   
import copy
import re
import warnings

def get_optimal_shuffle_encoding(folder_in, filename_in, encoding=None):
    """
    Function to find optimal shuffle encoding for a netcdf file.
    It will use deflate compression, and use trial and error to determine if shuffle should be used
    Encoding is an optional argument, to supply a scale_factor, dtype, _FillValue etc. NB: deflate compression and shuffle will be overwritten.
    """
    # Example
    # folder_in = r'O:\HybridDune experiment\data ADV, OBS\raw NetCDF\test'
    # filename_in = r'ADV_RWS1_Deployment1.nc'
    # encoding_shuffle_optimal = get_optimal_shuffle_encoding(folder_in, filename_in)
    # print(encoding_shuffle_optimal)

    ds = xr.open_dataset(os.path.join(folder_in, filename_in))
    filename_1var_shuffle_off = r'temp 1var shuffle_off' 
    filename_1var_shuffle_on  = r'temp 1var shuffle_on'  

    # Define two dictionaries for compression settings, one with shuffle on and one with shuffle off
    compression_shuffle_off = {var: {"zlib": True, "complevel": 4, 'shuffle': False} for var in list(ds.data_vars) + list(ds.coords)}  
    compression_shuffle_on  = {var: {"zlib": True, "complevel": 4, 'shuffle': True} for var in list(ds.data_vars) + list(ds.coords)}      

    # ---------------------------------------------------------------------------------------------------
    # Define two encoding dictionaries. Check if encoding is provided, otherwise use default compression 
    if encoding is None:
        encoding_shuffle_off = compression_shuffle_off
        encoding_shuffle_on  = compression_shuffle_on
    else:
        if list(encoding.keys())[0] != 'unlimited_dims':  # If useful encoding is found in the dataset
            warnings.warn( 
                f"Encoding starts with 'unlimited dims:', it is likely not useful (a default encoding created by XArray). In case of error, try not supplying encoding."
                f"Encoding: {encoding}"
            )

        # Combine supplied encoding with compression settings
        encoding_shuffle_off = copy.deepcopy(encoding)  # make a copy of the encoding dict, to add the shuffle setting
        encoding_shuffle_on  = copy.deepcopy(encoding)
        #print (f'Encoding supplied: {encoding}')  # print the encoding to check

        # Then add deflate compression to all variables and coordinates in netCDF, without overwriting existing keys
        for var, comp in compression_shuffle_off.items():  # for each variable in the dataset, 
            if var in encoding_shuffle_off:                # if the variable already has an encoding, update it with the compression settings
                encoding_shuffle_off[var].update(comp)
            else:                                          # if the variable does not have an encoding yet, add it 
                encoding_shuffle_off[var] = comp

        for var, comp in compression_shuffle_on.items():  # repeat for encoding_shuffle_on
            if var in encoding_shuffle_on:                
                encoding_shuffle_on[var].update(comp)
            else:                                          
                encoding_shuffle_on[var] = comp

    # ---------------------------------------------------------------------------------------------------
    # Now save the dataset with both shuffle settings, 1 variable at a time, compare the file sizes, then select optimal settings
    encoding_shuffle_optimal = copy.deepcopy(encoding_shuffle_on)  # make a copy of the encoding dict, to add the shuffle setting

    for var in ds.data_vars:        # for every variable in the dataset, if the number of elements in the variable is larger than 1000, determine optimal shuffle settings
        if ds[var].size > 1000: 
            #print(f'Processing variable: {var} with size {ds[var].size}')
            
            # make a new dataset, containing only the selected variable, dropping all other variables
            ds_1var = ds[[var]].copy()  # create a new dataset with only the selected variable

            # make a new encoding dictionary, with only variables and coordinates that are in ds_1var
            encoding_1var_shuffle_off = {var: encoding_shuffle_off[var] for var in list(ds_1var.data_vars) + list(ds_1var.coords)}
            encoding_1var_shuffle_on  = {var: encoding_shuffle_on[var]  for var in list(ds_1var.data_vars) + list(ds_1var.coords)}       
            
            # Note: some int vars give error when saving with scale_factor, in that case we need to cast them to float64 first
            try:
                ds_1var.to_netcdf(os.path.join(folder_in, filename_1var_shuffle_off), encoding=encoding_1var_shuffle_off)
            except Exception as e:
                if str(e).startswith("Cannot cast ufunc 'divide' output from dtype"):
                    print(f"Error for variable {var}: {e}. Casting variable to float64 and retrying.")
                    ds_1var[var] = ds_1var[var].astype('float64')
                    ds_1var.to_netcdf(os.path.join(folder_in, filename_1var_shuffle_off), encoding=encoding_1var_shuffle_off)
                    ds[var] = ds[var].astype('float64')  # also update the original dataset with the new data type. (order: after saving, so only if it works)
                else:
                    raise
                            
            ds_1var.to_netcdf(os.path.join(folder_in, filename_1var_shuffle_on),  encoding=encoding_1var_shuffle_on)

            # read the file size of the saved files
            size_shuffle_off = os.path.getsize(os.path.join(folder_in, filename_1var_shuffle_off))
            size_shuffle_on  = os.path.getsize(os.path.join(folder_in, filename_1var_shuffle_on))

            # if size_shuffle_off is smaller than size_shuffle_on, update the encoding_shuffle_optimal dictionary for the variable
            if size_shuffle_on < size_shuffle_off:
                print(f'Variable {var}: shuffle on is smaller: {size_shuffle_on/1024:.0f} kB < {size_shuffle_off/1024:.0f} kB')
            else:           
                encoding_shuffle_optimal[var] = encoding_shuffle_off[var]
                print(f'Variable {var}: shuffle on is larger: {size_shuffle_on/1024:.0f} kB > {size_shuffle_off/1024:.0f} kB')

    # Remove the temporary nc files
    os.remove(os.path.join(folder_in, filename_1var_shuffle_off))
    os.remove(os.path.join(folder_in, filename_1var_shuffle_on))

    return encoding_shuffle_optimal

