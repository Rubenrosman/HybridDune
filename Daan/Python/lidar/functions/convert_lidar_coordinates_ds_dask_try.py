# Dask-lazy version: convert lidar coordinates in chunks

# NB: Excessively slow. Use non-dask version, and load large datasets manually in parts if needed. 

import xarray as xr
import numpy as np

# Dask-aware sph2cart
def sph2cart(azimuth, elevation, r):
    rcoselev = r * np.cos(elevation)
    x = rcoselev * np.cos(azimuth)
    y = rcoselev * np.sin(azimuth)
    z = r * np.sin(elevation)
    return x, y, z

def convert_lidar_coordinates_ds_dask(ds, coordinate_system='local'):
    """
    Dask-lazy version: Convert lidar polar coordinates to Cartesian coordinates in either RD or local coordinate system.
    Keeps all large variables chunked/lazy for efficient memory use.
    """
    if coordinate_system not in ['local', 'RD']:
        raise ValueError("coordinate_system must be either 'local' or 'RD'")

    # Use xarray's astype for lazy conversion (dask-aware)
    profile_angle = ds.profile_angle.astype(np.float32)
    beam_angle = ds.beam_angle.astype(np.float32)
    r = ds.radius_lidar.astype(np.float32)  # This stays lazy/chunked

    # For small coordinate variables, .values is fine
    rotation_matrix = ds.rotation_matrix_lidar_to_RD.values.astype(np.float32)  # shape: (3, 3)
    translation_matrix = (
        np.array([ds.x_i_RD.values, ds.y_i_RD.values, ds.z_i.values]) - np.array([72400, 452100, 0])
    ).astype(np.float32)

    x_lidar, y_lidar, z_lidar = sph2cart(
        np.deg2rad(beam_angle), np.deg2rad(profile_angle), r
    )

    # Stack using xarray (keeps dask chunks)
    xyz_lidar = xr.concat([x_lidar, y_lidar, z_lidar], dim='cart').transpose(*r.dims, 'cart')
    xyz_lidar = xyz_lidar.chunk({'cart': -1})
    
    # Matrix multiplication: use xarray's apply_ufunc for dask compatibility
    xyz_RD = xr.apply_ufunc(
        lambda xyz: np.matmul(xyz, rotation_matrix) + translation_matrix,
        xyz_lidar,
        input_core_dims=[['cart']],
        output_core_dims=[['cart']],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[np.float32]
    )

    if coordinate_system == 'RD':
        ds_out = ds.copy(deep=True)
        ds_out = ds_out.drop_vars('radius_lidar')
        ds_out['x'] = (r.dims, (xyz_RD.sel(cart=0) + np.float32(72400)).data)
        ds_out['y'] = (r.dims, (xyz_RD.sel(cart=1) + np.float32(452100)).data)
        ds_out['z'] = (r.dims, xyz_RD.sel(cart=2).data)
        return ds_out

    else:  # coordinate_system == 'local'
        theta = np.deg2rad(np.float32(36.0))
        rotation_matrix_RDtoLocal = np.array([
            [np.cos(theta), np.sin(theta), 0],
            [-np.sin(theta), np.cos(theta), 0],
            [0, 0, 1]
        ], dtype=np.float32)
        shift = (np.array([71683.584, 452356.055, 0]) - np.array([72400, 452100, 0])).astype(np.float32)
        xyz_local = xr.apply_ufunc(
            lambda xyz: np.matmul(xyz - shift, rotation_matrix_RDtoLocal),
            xyz_RD,
            input_core_dims=[['cart']],
            output_core_dims=[['cart']],
            vectorize=True,
            dask='parallelized',
            output_dtypes=[np.float32]
        )
        ds_out = ds.copy(deep=True)
        ds_out = ds_out.drop_vars('radius_lidar')
        ds_out['x'] = (r.dims, xyz_local.sel(cart=0).data)
        ds_out['y'] = (r.dims, xyz_local.sel(cart=1).data)
        ds_out['z'] = (r.dims, xyz_local.sel(cart=2).data)
        return ds_out