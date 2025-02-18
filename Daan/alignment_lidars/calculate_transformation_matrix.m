clear, clc
xyz = [-1.3276	14.2143	-1.8844  % local coordinates of ground control points, storm 2, lidar 2, in m
       -1.3719	15.6844	1.3844   % file Lidar2__Pcl_0023_12_2024__12_06_30.062.PCD
       -1.3014	13.3857	2.952
       3.1823	-9.3287	-1.301
       3.0421	-7.2276	1.0061
       3.207	-9.6959	4.2199
       1.4108	8.0448	1.7928
       0	    0	    0];

xyz_RD_measured = [72488.47	    452098.15	5.45   % RD coordinates of ground control points
                   72487.617	452094.633	5.505
                   72484.851	452094.844	5.324
                   72469.414	452112.26	1.774
                   72469.66	    452109.102	1.841
                   72465.737	452108.107	1.706
                   72481.243	452099.093	2.879
                   72476.073	452105.429	4.692];
xyz_RD_measured = xyz_RD_measured - [72000 452000 0]  % subtract 72000 and 452000

% output cloudcompare
transformation_matrix_CC = [-0.019893104	 0.789837182	-0.612993896	476.05029
                             0.032912239	-0.612265587	-0.789966881	105.4001
                            -0.999260247	-0.035889894	-0.013815446	4.627503872];

rotation_matrix_CC = transformation_matrix_CC(1:3,1:3)
translation_matrix_CC = transformation_matrix_CC(1:3,4)

% calculate transformation matrix directly from coordinates
[rotation_matrix_calc, translation_matrix_calc, lrms] = Kabsch(xyz', xyz_RD_measured')  % 3x3 rotation matrix; 3x1 translation matrix; RMSE. Equal to output cloudcompare.

% check quality
rotation_matrix_CC - rotation_matrix_calc;        % same result as cloud compare
translation_matrix_CC - translation_matrix_calc;  % same result

xyz_new_calc = xyz * rotation_matrix_calc' + translation_matrix_calc'
error = xyz_RD_measured - xyz_new_calc; % difference measured and calculated RD x,y,z coordinates, in m
distance = sqrt( error(:,1).^2 + error(:,2).^2 + error(:,3).^2 ); % total error per point, in m