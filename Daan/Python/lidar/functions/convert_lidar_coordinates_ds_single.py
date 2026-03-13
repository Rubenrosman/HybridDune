# Import necessary libraries --------------------------------------------------------------
import numpy as np
import xarray as xr
import pandas as pd
from warnings import warn

# define helper functions ------------------------------------------------------------------
def sph2cart(azimuth,elevation,r):
    rcoselev = r * np.cos(elevation)
    x = rcoselev * np.cos(azimuth)
    y = rcoselev * np.sin(azimuth)
    z = r * np.sin(elevation)
    return x, y, z

# Main function: convert dataset ----------------------------------------------------------------------
def convert_lidar_coordinates_ds(ds, coordinate_system='local', max_echos=False, keep_I=True, keep_r=True, keep_y=True):
    """
    Convert lidar polar coordinates to Cartesian coordinates in either RD or local coordinate system.
    Calculation and output in single precision (float32) for better performance.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Dataset containing lidar data in polar coordinates with variables:
        - profile_angle, beam_angle, radius_lidar
        - rotation_matrix_lidar_to_RD, x_i_RD, y_i_RD, z_i
    coordinate_system : str, optional
        Choice of output coordinate system: 'local' (default) or 'RD'
    max_echos : bool, optional
        If True, only the maximum radius and corresponding intensity along the 'echos' dimension are returned. Default is False.
    keep_I : bool, optional
        If True (default), intensity variable is kept in output dataset. 
    keep_r : bool, optional
        If True (default), radius_lidar variable is kept in output dataset. 
    keep_y : bool, optional
        If True (default), y variable is kept in output dataset. If false and coordinate_system is 'local', alongshore coordinate y is dropped.
        
    Returns:
    --------
    xarray.Dataset
        Dataset with Cartesian coordinates (x, y, z) replacing/in addition to radius_lidar

    Example usage:
    ---------------
    nc_in = r"O:\HybridDune experiment\data lidar netcdf\storm1_lidar1_polar_10sinterval - new.nc"

    t0    = pd.Timestamp('2024-12-18 14:10')
    t_end = pd.Timestamp('2024-12-18 18:01')
    ds = xr.open_dataset(nc_in).sel(t=slice(t0, t_end))   # select specific time range: t0 until t_end
    
    ds_local = convert_lidar_coordinates_ds(ds, 'local')
    ds_RD = convert_lidar_coordinates_ds(ds, 'RD')
    """
    
    # Input validation
    if coordinate_system not in ['local', 'RD']:
        raise ValueError("coordinate_system must be either 'local' or 'RD'")
    
    # Extract variables from dataset
    profile_angle      = ds.profile_angle.astype(np.float32)
    beam_angle         = ds.beam_angle.astype(np.float32)
    r                  = ds.radius_lidar.astype(np.float32)
    rotation_matrix    = ds.rotation_matrix_lidar_to_RD.values.astype(np.float32)  # shape: (3, 3)
    translation_matrix = (  np.array([ds.x_i_RD.values, ds.y_i_RD.values, ds.z_i.values])  -  np.array([72400, 452100, 0])  ).astype(np.float32) # shift coordinates: smaller numbers for better numerical precision
    
    if not keep_r:
        ds = ds.drop_vars('radius_lidar')  # NB: remove line if r ever needs to be returned in output

    if keep_I:
        I = ds.intensity.astype(np.float16) # Convert to float16, integers with max 255 so float32 not needed. Don't use int8, data contains nans.
    else:
        ds = ds.drop_vars('intensity', errors='ignore')

    # calculate maximum of r and I along 'echos' dimension, if asked for. 
    if max_echos:
        r = r.max(dim='echos')                                        # r is sorted over echos, so the max r is the last non-nan r
        if keep_I:
            ds['intensity']    = I.ffill(dim='echos').isel(echos=-1)  # take last non-nan intensity (= I of max r, = I of last r)
    else:
        if keep_I:
            ds['intensity'] = I

    if keep_r:
        ds['radius_lidar'] = r

    # Convert polar coordinates to cartesian coordinates, in local lidar coord system
    [x_lidar, y_lidar, z_lidar] = sph2cart(np.deg2rad(beam_angle), np.deg2rad(profile_angle), r)
    
    # Convert to RD system
    # Stack lidar coordinates into shape (..., 3) for matrix multiplication
    xyz_lidar = np.stack([x_lidar, y_lidar, z_lidar], axis=-1)
    
    # Apply rotation (matrix multiply last axis) and translation
    xyz_RD = np.matmul(xyz_lidar, rotation_matrix) + translation_matrix
    
    # Prepare output
    # ds = ds.drop_vars('radius_lidar')  # Moved up

    # if ds[y_i_local] does not exist yet (missing in current NC files), compute it
    if 'y_i_local' not in ds:
        # add local y-coordinate to dataset (missing)
        x_i_RD = ds['x_i_RD'].values
        y_i_RD = ds['y_i_RD'].values

        xy_RD = np.array([x_i_RD, y_i_RD]).T
        theta = np.deg2rad(36)
        transformation_matrix = np.array([ [np.cos(theta), np.sin(theta)],[-np.sin(theta), np.cos(theta)] ])
        xy_loc = ( xy_RD - [71683.584, 452356.055] ) @ transformation_matrix
        x_i_loc = xy_loc.T[0]
        y_i_loc = xy_loc.T[1]

        ds['y_i_local'] = y_i_loc

    if coordinate_system == 'RD':
        # Assign RD coordinates to output ds
        ds['x'] = (r.dims, xyz_RD[..., 0] + np.float32(72400))   # +72400: correct the shift
        ds['y'] = (r.dims, xyz_RD[..., 1] + np.float32(452100))
        ds['z'] = (r.dims, xyz_RD[..., 2])
        if not keep_y:
            warn("keep_y is ignored when coordinate_system is 'RD', y-coordinate is always kept.")
    
    else:  # coordinate_system == 'local'
        # Transform to x,y coordinates in local coordinate system (alongshore, cross shore)
        theta = np.deg2rad(np.float32(36.0))  # local grid is rotated 36 degrees clockwise wrt RD system
        rotation_matrix_RDtoLocal = np.array([
            [np.cos(theta), np.sin(theta), 0],
            [-np.sin(theta), np.cos(theta), 0],
            [0, 0, 1]
        ], dtype=np.float32)
        
        # Subtract zero-point of local grid, then rotate 36 degrees
        # xyz_local = xyz_RD - np.array([71683.584, 452356.055, 0], dtype=np.float32)  # normal version
        xyz_local = xyz_RD - (  np.array([71683.584, 452356.055, 0])  -  np.array([72400, 452100, 0])  ).astype(np.float32) # for shifted RD coordinates, better numberical precision
        xyz_local = np.matmul(xyz_local, rotation_matrix_RDtoLocal)
        
        # Assign local coordinates to output ds
        ds['x'] = (r.dims, xyz_local[..., 0])
        if keep_y:
            ds['y'] = (r.dims, xyz_local[..., 1])
        ds['z'] = (r.dims, xyz_local[..., 2])
        
    return ds