%% plot 4 specific profiles of two different storms
clear,clc
figure(2), clf(2), hold on

for n=1%:4
map{1} = 'O:\2024-12-18 to 2024-12-20, Storm 1\Lidars\20241220_LiDAR2\10s interval data\';  % storm 1
map{2} = 'O:\2024-12-18 to 2024-12-20, Storm 1\Lidars\20241220_LiDAR2\10s interval data\';  % storm 1
map{3} = 'O:\2024-12-18 to 2024-12-20, Storm 1\Lidars\20241220_LiDAR2\10s interval data\';  % storm 1
map{4} = 'O:\2025-01-06 to 2025-01-07, Storm 5\Lidars\20250108_LiDAR2\2025_01_06_13hr\';    % storm 5

naam{1} = 'Lidar2__Pcl_0018_12_2024__14_25_30.214.PCD';   % t1: 18 dec, 14:25 (low tide, just before storm
naam{2} = 'Lidar2__Pcl_0019_12_2024__13_07_20.214.PCD';   % t2: day later, again low tide
naam{3} = 'Lidar2__Pcl_0019_12_2024__18_07_20.214.PCD';   % t3: half day later, high tide
naam{4} = 'lidar2__Pcl_0006_01_2025__13_00_00.024.PCD';   % t4: storm 2

map_naam = [map{n},naam{n}];
ptCloud = pcread(map_naam); % load pointcloud 

% convert to polar coordinates
xyz = ptCloud.Location;
intensity = ptCloud.Intensity;
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
% [azimuth,elevation,r] = cart2sph(x,y,z);
% 
% % set polar coordinates in degrees and sort points
% azimuth  = rad2deg(azimuth );
% elevation = rad2deg(elevation);
% 
% % sort rows by profile number, then angle coordinate
% elevation_round = round(elevation,2);
% sphere_coor = [elevation_round,azimuth,elevation,r,intensity]; % temporary add rounded elevation coord. for sorting
% sphere_coor = sortrows(sphere_coor,[1 2]);  % sort first by the first column (rounded elev.), then the second (azi.)
% 
% % determine profile number 
% n_profile = single( [1; 1 + cumsum(logical(diff(sphere_coor(:,1))))] );
% 
% % detemrine angle number
% n_hoek = ceil( (sphere_coor(:,2)+180)*2 );
% 
% % determine echo count
% azimuth_round = round( sphere_coor(:,2), 2 );
% echo_count = count_consequetive_equal_down(azimuth_round);
% 
% sphere_coor = sphere_coor(:,2:end); % remove first column (rounded elevation coordinates)
% 
% % save as new pointcloud variable
% azimuth  = deg2rad(sphere_coor(:,1));
% elevation = deg2rad(sphere_coor(:,2));
% 
% [x,y,z] = sph2cart(azimuth,elevation,sphere_coor(:,3));
% xyz = [x,y,z];
% I   = sphere_coor(:,4);
% 
% i_recht = n_profile==3;
x2 = x/1000;
y2 = y/1000;
z2 = z/1000;
xyz2=[x2,y2,z2];


scatter(y2,z2)
end
legend

%% Load .mat file with 10s data of full storm 

clear,clc
load('O:\2024-12-18 to 2024-12-20, Storm 1\Lidars\20241220_LiDAR2\10s interval data.mat')  % load file with 10s data of storm1, lidar2

%% plot from above dataset
figure(1), clf(1), hold on
t=5000:72:6000;   % select time steps to use. here every 72 files, so every 72*10s=720 s. i.e. 5 plots per hour
color_map = colorcet('GOULDIAN','N',length(t));  % colormap, to assign different times different colors 
% color_map = parula(length(t)); % colorcet is downloaded function. For native matlab alternative, comment line above, uncomment this one

for n = 1:length(t)        % for selected time steps
    for profiel = 3%1:16   % and selected profiles.  profile 3 is perdendicular to coastline
        xyz = xyz_raw_mm{t(n),profiel}; % select xyz
        x = xyz(:,2) / 1000;  % cross shore. convert to m
        y = xyz(:,3) / 1000;  % alongshore
        z = -xyz(:,1) / 1000; % vertical

        scatter(x,z,'.','MarkerEdgeColor',color_map(n,:))
    end
    labels{n} = string(datetime_all(t(n)),'dd-MM-yy HH:mm');  % make label from date-time of point cloud
end
xlim([-5 15])
ylim([-4 1])
legend(labels) % make legend with date-time labels
%% fucntion
function n_con = count_consequetive_equal_down(A)
% Count consequetive ones in matrix. Count from top to bottom. Output of int16 class. Works for vectors and matrices. Example:
% 1 0 1 1  gives   1 0 1 1
% 0 1 1 1          0 1 2 2
M = size(A,1);

% where does A change, i.e., where does a new run start
Q = A ~= circshift(A,[1 0]); 
Q(1,:,:) = 1;  % first row: change by definition

Q_index     = Q .* ((0:M-1)'); % which column did it change
Q_index_max = cummax(Q_index,1);     % what is the last column that changed?

n_con       = (1:(M))' - Q_index_max; % How long is the series? So current column - last column
end