clear, clc
% V2: save as polar coordinates
% V3: save as 3d matrix
% V3b: save full folder. Comment code that is not used. 
% V4: save as video
% V5: as V3b, cleaned
map_all      = {'O:\2024-12-22, Storm 2\Lidars\20241223_LiDAR2\10s interval data\';
                'O:\2024-12-22, Storm 2\Lidars\20241223_LiDAR3\10s interval data\';
                'O:\2025-01-01, Storm 3\Lidars\20250102_LiDAR2\10s interval data\';
                'O:\2025-01-01, Storm 3\Lidars\20250102_LiDAR3\10s interval data\';
                'O:\2025-01-06 to 2025-01-07, Storm 5\Lidars\20250108_LiDAR5\10s interval data\'};

map_save_all = {'O:\2024-12-22, Storm 2\Lidars\20241223_LiDAR2\';
                'O:\2024-12-22, Storm 2\Lidars\20241223_LiDAR3\';
                'O:\2025-01-01, Storm 3\Lidars\20250102_LiDAR2\';
                'O:\2025-01-01, Storm 3\Lidars\20250102_LiDAR3\';
                'O:\2025-01-06 to 2025-01-07, Storm 5\Lidars\20250108_LiDAR5\'};

for n_map = 1:2
    clearvars -except map_all map_save_all n_map
    disp(n_map)

    map = map_all{n_map};
    map_save = map_save_all{n_map};

%% maak lijst met bestandsnamen
%     map = mappen_in{n_map}
%     mkdir([map,'10s interval data\'])
%     mkdir([mappen_uit{n_map},'10s interval data\'])
list_files = dir(map);
disp('step 1 finished')

%% 2 Bepaal tijd uit bestand
n_1=1;
while list_files(n_1).isdir || ~strcmpi( list_files(n_1).name(end-2:end), 'PCD') % if entry is a folder, or wrong file type, 
    n_1=n_1+1;  % try next file
end

n_max = length(list_files);
while list_files(n_max).isdir || ~strcmpi( list_files(n_max).name(end-2:end), 'PCD') % if entry is a folder, or wrong file type, 
    n_maxn_max - 1;  % try next file
end

% n_max = 1000;             % optional, for testing code. 
% warning('n_max is 1000')

naam1 = list_files(n_1).name;   % bepaal eerste bestandsnaam. Eerst 'lege bestanden'.

% find datum uit bestand
n_letter_date = strfind(naam1, 'Pcl') + 6;
date_string1 = naam1(n_letter_date : n_letter_date + 23);
t1 = datetime(date_string1, 'InputFormat','dd_MM_yyyy__HH_mm_ss.SSS');

% make lijst met alle data
datetime_all = repmat(t1, n_max - n_1+ 1,1);
parfor n = n_1 : n_max
    %     naam = lijst_bestanden(n).name;
    if length(list_files(n).name) ~= length(naam1)
        error('different file format: file name has different length')
    end
    datetime_all(n,1) = datetime(list_files(n).name(n_letter_date : n_letter_date + 23), 'InputFormat','dd_MM_yyyy__HH_mm_ss.SSS');
end
disp('step 2 finished')


%% Load pointclouds
tic
parfor n = n_1 : n_max
    map_naam = [map,list_files(n).name];
    [r_reshaped, I_reshaped,hoek_profile,hoek_ini_at_profile, xyz_cell,I_cell] = load_pointcloud_2(map_naam);
    
    xyz_cell_all(n,:) = xyz_cell;
    I_cell_all(n,:)   = I_cell;

    angle_profile_all(n,:) = hoek_profile;
    angle_ini_at_profile_all(n,:) = hoek_ini_at_profile;
end
disp('step 3 finished')
toc

%% skip non-existing n
datetime_all = datetime_all(n_1:n_max);
list_files = list_files(n_1:end);

angle_profile_all = angle_profile_all(n_1:n_max,:);
angle_ini_at_profile_all = angle_ini_at_profile_all(n_1:n_max,:);

xyz_cell_all = xyz_cell_all(n_1:n_max,:);
I_cell_all = I_cell_all(n_1:n_max,:);

% simplify
if max(abs(diff(angle_profile_all))) < 0.0001
    angle_profile = angle_profile_all(1,:);
end
if max(abs(diff(angle_ini_at_profile_all))) < 0.0001
    angle_ini_at_profile = angle_ini_at_profile_all(1,:);
end

%% save


% rename outputs
xyz_raw_mm = xyz_cell_all;
intensity_256 = I_cell_all;

save([map_save,'10s interval data.mat'], 'xyz_raw_mm','intensity_256', 'datetime_all','list_files','-v7.3')
% save([map_save,'10s interval data.mat'], 'xyz_raw_mm','intensity_256', 'datetime_all','list_files','-v7.3') % only use if variable is too large for normal ('-v7'): 10%-20% larger files

disp('saved')

end
%% fucntion load_pointcloud_1.      not used in this script. not checked if it is still correct.
function [x,y,z,intensity,n_profile] = load_pointcloud(map_naam)
% Load point cloud ----------------------------------------------------------------------------------------------------
ptCloud = pcread(map_naam); % load pointcloud

x = ptCloud.Location(:,1);
y = ptCloud.Location(:,2);
z = ptCloud.Location(:,3);
intensity = ptCloud.Intensity;

[angle_laserbeam,angle_laserplane,r] = cart2sph(x,y,z);

% set polar coordinates in degrees and sort points
angle_laserbeam  = rad2deg(angle_laserbeam );
angle_laserplane = rad2deg(angle_laserplane);

% Define points by profile number and obs_number within profile -------------------------------------------------------
% sort rows by profile number, then angle coordinate
sphere_coor = [ angle_laserplane, angle_laserbeam, r, intensity]; 
[sphere_coor_sorted, sort_index] = sortrows_rounded(sphere_coor, 2);  % sort rows (first by laser plane angle, then by laser beam angle). ', 2':  
[~,sort_index_inverse] = sort(sort_index);  % used later to reverse sorting and check results

% determine profile number.     n_profile is a nubmer between 1 and 16 (16 profiles)
angle_laserplane_sorted_rounded = round(sphere_coor_sorted(:,1), 2); % round angles of laser plane. Ensures that changes in angle are really a different plane, not just numerial accuracy
n_profile = [1; 1 + cumsum(logical(diff(angle_laserplane_sorted_rounded)))]; % logical(diff(coordinate): if the (rounded) coordinate changes, a new profile starts. Cumsum: count changes
if max(n_profile) > 16,     error('calcualted max profile number is more than 16'),     end

% determine angle number, i.e. observation no. within profile
% NB: the Sick Multiscan lidar measures exactly every 0.5 degrees. So instead of recording the exact angle of a
% measurement, we can simply record the that it is the nth measurement within that profile.
n_angle_within_profile = ceil( (sphere_coor_sorted(:,2)+180)*2 ); % number between 1 and 720. An observation between 0°-0.5° gives n_angle=1. 0.5°-1° gives 2. etc.


% save as new pointcloud variable
angle_laserplane = deg2rad(sphere_coor(:,1));
angle_laserbeam  = deg2rad(sphere_coor(:,2));

[x,y,z] = sph2cart(angle_laserbeam,angle_laserplane,sphere_coor(:,3));
end

% function 2: polar coordinates ---------------------------------------------------------------------------------------
%% fucntion load_pointcloud_2
function [r_reshaped,  I_reshaped,angle_profile,angle_ini_at_profile, xyz_cell,I_cell] = load_pointcloud_2(map_naam)
% Load point cloud ----------------------------------------------------------------------------------------------------
ptCloud = pcread(map_naam); % load pointcloud

x = ptCloud.Location(:,1);
y = ptCloud.Location(:,2);
z = ptCloud.Location(:,3);
intensity = ptCloud.Intensity;

[angle_laserbeam,angle_laserplane,r] = cart2sph(x,y,z);

% set polar coordinates in degrees and sort points
angle_laserbeam  = rad2deg(angle_laserbeam );
angle_laserplane = rad2deg(angle_laserplane);

% Define points by profile number and obs_number within profile -------------------------------------------------------
% sort rows by profile number, then angle coordinate
sphere_coor = [ angle_laserplane, angle_laserbeam, r, intensity]; 
[sphere_coor_sorted, sort_index] = sortrows_rounded(sphere_coor, 2);  % sort rows (first by laser plane angle, then by laser beam angle). ', 2':  
[~,sort_index_inverse] = sort(sort_index);  % used later to reverse sorting and check results

% determine profile number.     n_profile is a nubmer between 1 and 16 (16 profiles)
angle_laserplane_sorted_rounded = round(sphere_coor_sorted(:,1), 2); % round angles of laser plane. Ensures that changes in angle are really a different plane, not just numerial accuracy
n_profile = [1; 1 + cumsum(logical(diff(angle_laserplane_sorted_rounded)))]; % logical(diff(coordinate): if the (rounded) coordinate changes, a new profile starts. Cumsum: count changes
if max(n_profile) > 16,     error('calcualted max profile number is more than 16'),     end

% determine angle number, i.e. observation no. within profile
% NB: the Sick Multiscan lidar measures exactly every 0.5 degrees. So instead of recording the exact angle of a
% measurement, we can simply record the that it is the nth measurement within that profile.
n_angle_within_profile = ceil( (sphere_coor_sorted(:,2)+180)*2 ); % number between 1 and 720. An observation between 0°-0.5° gives n_angle=1. 0.5°-1° gives 2. etc.

% determine echo count. Count no. of obs at the same laser angle. (echos could be reflections, or maybe more than a full laser rotation)
echo_count = count_consequetive_equal_down_margin(sphere_coor_sorted(:,2), 0.01); % function to count how many obs in consequetive rows have the same value. same within 0.01 deg, to account for numerical inaccuracy

% define variables to save
I_int8 = uint8(sphere_coor_sorted(:,4) * 255);  % store intensity as integer value, between 0 and 255. Requires less memory (but less precize)
r_int16 = uint16(sphere_coor_sorted(:,3));      % store r as integer value. can handle values up to 65535 mm, so 65 m. (Sick Multiscan 165 has 62 m range, interally 1mm accuracy for range, so no cutting nor rounding should occur)
if max(r) > uint16(inf),    error('Integer overflow. Max distance is more than 65535, cannot be stored as uint16 class'),   end

% put r, I in full matrix (720 values x 16 beams x 3 echos
r_reshaped = zeros(720,16,3,'uint16');  % initialize as zeros: where there is no measurement, distance r will be zero
I_reshaped = zeros(720,16,3,'uint8');

for n=1:length(r)
    r_reshaped( n_angle_within_profile(n),n_profile(n),echo_count(n) ) = r_int16(n);
    I_reshaped( n_angle_within_profile(n),n_profile(n),echo_count(n) ) = I_int8(n);
end

% output r and I as cell
xyz_out = [x(sort_index),y(sort_index),z(sort_index)];
for n=1:16
    xyz_cell{1,n} = xyz_out(n_profile == n,:);
    I_cell{1,n} = I_int8(n_profile == n);
end

%     % echos los opslaan
%     a = [sphere_coor_sorted,echo_count];
%     sphere_coor_echos=[];
%     for n=1:length(a)
%         if echo_count(n)>1
%             sphere_coor_echos = [sphere_coor_echos; n_hoek(n) n_profile(n) r_sorted(n) I(n)*255 echo_count(n)];
%         end
%     end
%     sphere_coor_echos_uint16 = uint16(sphere_coor_echos);
    
% save exact angles
n_profile_start = [1; find(diff(n_profile))+1]; % find at which indices a new profile starts
angle_profile = sphere_coor_sorted(n_profile_start,1)';
angle_ini_at_profile = [sphere_coor_sorted(n_profile_start,2) - ( n_angle_within_profile(n_profile_start) - 1 ) * 0.5]';

% Check results: is compressed point cloud identical to original? ------------------------------------------------------
% Calculate polar coordinates from distance matrix 'r', and comvert to xyz coordinates
r_check = single(r_reshaped);  % everything single: allow for floating point numbers
angle_laserplane_check = single( repmat( angle_profile, 720, 1, 3) );
angle_laserbeam_check  = single( repmat(angle_ini_at_profile,720,1,3) + [0:0.5:359.9]' );

[x_check,y_check,z_check] = sph2cart( deg2rad(angle_laserbeam_check), deg2rad(angle_laserplane_check), r_check ); % convert to xyz coordinates, in mm
xyz_check = cat(4,x_check,y_check,z_check);                            % put together in single matrix

% Put xyz coordinates in same order as in original pointcloud file. Check: identical?
r_backshaped = reshape( permute(r_reshaped,[3 1 2]),   720*3*16,1);     % reshape to a single column vector, similar to vector r.   permute: ensure correct order ...
xyz_check    = reshape( permute(xyz_check,[3 1 2 4]), 720*3*16,3);     %    of elements. first echo 1-3 of profile_1, angle_1. then echo 1-3 of angle_2. etc...

xyz_check    = xyz_check( find(r_backshaped), : );                     % remove all the rows where r is zero (no measurement)
xyz_check    = xyz_check(sort_index_inverse,:);                        % revert sorting done at start, so that xyz can be compared to original pointcloud
xyz_ori = ptCloud.Location;                                            % xyz from original pointcloud
if max(abs(xyz_check-xyz_ori),[],'all') > 0.1
    max(abs(xyz_check-xyz_ori))
    error('The difference between the original and compressed pointcloud is more than 0.1 mm')
end
end

%% function consequetive equals -----------------------------------------------------------------------------------------
function n_con = count_consequetive_equal_down_margin(A,margin)
% Count consequetive equal values in matrix. Count from top to bottom. Works for vectors and matrices. Equal within
% margin

equal_within_margin = [zeros(1, size(A,2) );  abs( A(1:end-1,:) - A(2:end,:)) < margin];
M = size(equal_within_margin,1);

% where does A change
Q = equal_within_margin ~= circshift(equal_within_margin,[1 0]); 
Q(1,:,:) = 1;  % first row: change by definition

Q_index     = Q .* ((0:M-1)'); % if it changes (Q=1), which row did it change
Q_index_max = cummax(Q_index,1);     % what is the last row that changed?

n_con       = (1:(M))' - Q_index_max; % How long is the series? So current column - last column
n_con       = 1 + n_con .* equal_within_margin;     % Only count series of trues, make 0 for false
end

%% function sortrows_rounded
function [matrix_sorted, sort_index] = sortrows_rounded (matrix_in,decimals)
% Function to sort matrix, while rounding almost equal values to make them identical
% e.g. rounding to 0.1 (decimals=1) ensures that a matrix of A = [1.10 2 1
%                                                                 1.14 3 2
%                                                                 1.16 1 3
%                                                                 2    1 4]
% is sorted in row order 3-1-2-4, whereas sortrows(round(A,1)) whould have given 2-1-3-4

tolerance = 10^-decimals;

% sort, round first column where closer together than tolerance
[matrix_sorted, sort_index] = sortrows(matrix_in);          % sorting
matrix_round_sorted = matrix_sorted;

for n = 2 : size(matrix_in,1)                               % rounding first column where closer together than tolerance
    if abs( matrix_round_sorted(n,1)-matrix_round_sorted(n-1,1))< tolerance
        matrix_round_sorted(n,1) = matrix_round_sorted(n-1,1);
    end
end
[~,sort_index_inverse] = sort(sort_index);                  % make rounded matrix, but in order of original matrix
matrix_round = matrix_round_sorted(sort_index_inverse,:);

% repeat with second column: sort, round second column where diff < tolerance. 
[matrix_sorted, sort_index] = sortrows(matrix_round);       % sorting again. NB: obtained rounded matrix, so that sorting 
matrix_round_sorted = matrix_sorted;                        %      over first column is 'correct'

for n = 2 : size(matrix_in,1)                               % rounding second column where diff < tolerance
    if abs( matrix_round_sorted(n,2)-matrix_round_sorted(n-1,2))< tolerance
        matrix_round_sorted(n,2) = matrix_round_sorted(n-1,2);
    end
end
[~,sort_index_inverse] = sort(sort_index);                  % again: make rounded matrix, but in order of original matrix
matrix_round = matrix_round_sorted(sort_index_inverse,:);

% Obtain final sorting order, sorted matrix
[~, sort_index] = sortrows(matrix_round); % use rounded matrix for correct sorting order. 
matrix_sorted = matrix_in(sort_index,:);  % apply order to original (un-rounded) matrix and output result
end