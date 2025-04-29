function [xyz_local, xyz_RD] = convert_lidar_coordinates(xyz_lidar,n_lidar, n_storm)
% Lidar transformation for hybrid dune project. From xyz in lidar coordinates (in m), to local project coordinates or RD coordinates
% xyz in meters, matrix with 3 columns (n rows, n highier dimensions) 

% Load transformation matrix for conversion to RDNAP
% load from excel
sheetname = ['storm',num2str(n_storm),' lidar',num2str(n_lidar)];
transformation_matrix = readmatrix('C:\Matlab\hybrid_dunes\lidar_transformation_matrices.xlsx','Sheet',sheetname,'Range','A1:D4');
% transformation_matrix = readmatrix('O:\HybridDune experiment\lidar_transformation_matrices.xlsx','Sheet',sheetname,'Range','A1:D4');

% split in rotation and translation, apply
rotation_matrix =  transformation_matrix(1:3,1:3)';
translation_matrix = transformation_matrix(1:3,4)' + [72000 452000 0];

% xyz_RD = xyz_lidar * rotation_matrix + translation_matrix;
xyz_RD = pagemtimes(xyz_lidar, rotation_matrix) + translation_matrix;

% Convert to local cross_shore, alongshore coordinates 
rotation_matrix2 = [ cosd(36) sind(36) 0
                    -sind(36) cosd(36) 0
                            0       0  1];

x0 = 71683.584;
y0 = 452356.055;
% xyz_local = (xyz_RD - [x0 y0 0]) * rotation_matrix2;
xyz_local = pagemtimes(xyz_RD - [x0 y0 0], rotation_matrix2);

end