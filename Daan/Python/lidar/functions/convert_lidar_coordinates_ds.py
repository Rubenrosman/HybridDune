# Import necessary libraries --------------------------------------------------------------
import numpy as np
import xarray as xr
import pandas as pd

# define helper functions ------------------------------------------------------------------
def sph2cart(azimuth,elevation,r):
    rcoselev = r * np.cos(elevation)
    x = rcoselev * np.cos(azimuth)
    y = rcoselev * np.sin(azimuth)
    z = r * np.sin(elevation)
    return x, y, z

# def rad2deg(angleInRadians):
#     angleInDegrees = 180/np.pi * angleInRadians
#     return angleInDegrees

# def deg2rad(angleInDegrees):
#     angleInRadians = np.pi/180 * angleInDegrees
#     return angleInRadians

# Main function: convert dataset ----------------------------------------------------------------------
def convert_lidar_coordinates_ds(ds, coordinate_system='local'):
    """
    Convert lidar polar coordinates to Cartesian coordinates in either RD or local coordinate system.
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Dataset containing lidar data in polar coordinates with variables:
        - profile_angle, beam_angle, radius_lidar
        - rotation_matrix_lidar_to_RD, x_i_RD, y_i_RD, z_i
    coordinate_system : str, optional
        Choice of output coordinate system: 'local' (default) or 'RD'
        
    Returns:
    --------
    xarray.Dataset
        Dataset with Cartesian coordinates (x, y, z) replacing radius_lidar

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
    translation_matrix = np.array([ds.x_i_RD.values, ds.y_i_RD.values, ds.z_i.values], dtype=np.float32)  # shape: (3,)
    
    # Convert polar coordinates to cartesian coordinates, in local lidar coord system
    [x_lidar, y_lidar, z_lidar] = sph2cart(np.deg2rad(beam_angle), np.deg2rad(profile_angle), r)
    
    # Convert to RD system
    # Stack lidar coordinates into shape (..., 3) for matrix multiplication
    xyz_lidar = np.stack([x_lidar, y_lidar, z_lidar], axis=-1)
    
    # Apply rotation (matrix multiply last axis) and translation
    xyz_RD = np.matmul(xyz_lidar, rotation_matrix) + translation_matrix
    
    if coordinate_system == 'RD':
        # Create xarray dataset with RD coordinates
        ds_out = ds.copy(deep=True)
        ds_out = ds_out.drop_vars('radius_lidar')
        ds_out['x'] = (r.dims, xyz_RD[..., 0])
        ds_out['y'] = (r.dims, xyz_RD[..., 1])
        ds_out['z'] = (r.dims, xyz_RD[..., 2])
        return ds_out
    
    else:  # coordinate_system == 'local'
        # Transform to x,y coordinates in local coordinate system (alongshore, cross shore)
        theta = np.deg2rad(np.float32(36.0))  # local grid is rotated 36 degrees clockwise wrt RD system
        rotation_matrix_RDtoLocal = np.array([
            [np.cos(theta), np.sin(theta), 0],
            [-np.sin(theta), np.cos(theta), 0],
            [0, 0, 1]
        ], dtype=np.float32)
        
        # Subtract zero-point of local grid, then rotate 36 degrees
        xyz_local = xyz_RD - np.array([71683.584, 452356.055, 0], dtype=np.float32)
        xyz_local = np.matmul(xyz_local, rotation_matrix_RDtoLocal)
        
        # Create xarray dataset with local coordinates
        ds_out = ds.copy(deep=True)
        ds_out = ds_out.drop_vars('radius_lidar')
        ds_out['x'] = (r.dims, xyz_local[..., 0])
        ds_out['y'] = (r.dims, xyz_local[..., 1])
        ds_out['z'] = (r.dims, xyz_local[..., 2])
        return ds_out