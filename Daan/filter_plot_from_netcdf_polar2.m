% Try plotting profiles
%% load data
clear, clc
map_file = 'O:\HybridDune experiment\2024-12-18 to 2024-12-20, Storm 1\Lidars\20241220_LiDAR1\storm1_lidar1_polar.nc';
% lidar/storm, also change lines for filtering (x_skip, y_skip) and which profiles to plot (n_profile_plot)
netcdf_info = ncinfo(map_file);

figure(4), clf(4), hold on
colorcet('R1') % set colors of plot, to ensure the colorbar hasthe same colorscale as the lines
lijn_kleuren = colorcet('R1','N',84); % define colors of the lines. Set a rainbow color scale with as many colors as no. of lines

for n_plot = 1:83
    n_start  = [1+n_plot*1800*4, 1,   3, 1];    % first pointcloud, first point in profile, third profile
    n_count  = [1800*4, 720, 1, 3];  % download over how many obs/gaps
    n_stride = [1,   1,   1, 1];  % with which gap width

    radius_lidar   = single( ncread(map_file,'radius_lidar',n_start,n_count,n_stride) )/1000; % radius, uint16
    profile_angle  = single( ncread(map_file,'profile_angle', n_start(3),1) )'; % profile angle, 1x16
    ini_beam_angle = single( ncread(map_file,'ini_beam_angle',n_start(3),1) )'; % profile angle, 1x16
    t_string       =   ncread(map_file,'time_string',n_start(1),n_count(1),n_stride(1));

    angle_laserplane = reshape(profile_angle,1,1,[]);  % in third dimension: third dim om radius is the proifles
    angle_laserbeam  = reshape(ini_beam_angle,1,1,[]); % idem
%     angle_laserplane = repmat(angle_laserplane, n_count(1), 720, 1);  % repeat: make T x 720 x n_profiles x 3
    angle_laserbeam  = repmat(angle_laserbeam, 1, 720, 1) + [0:0.5:359.9];

    %% Check if results are the same
    %{
[x_check,y_check,z_check] = sph2cart( deg2rad(angle_laserbeam), deg2rad(angle_laserplane), r ); % convert to xyz coordinates, in mm
xyz_check = cat(5,x_check,y_check,z_check);                            % put together in single matrix

% Put xyz coordinates in same order as in original pointcloud file. Check: identical?
r_backshaped = reshape( permute(radius_lidar,[4 2 3 1]),   720*3*n_count(3),1,n_count(1));  % reshape to a single column vector per pointcloud (720x3x"16" column) x 1 x T.  permute: ensure correct order ...
xyz_check    = reshape( permute(xyz_check,[4 2 3 5 1]), 720*3*n_count(3),3,n_count(1));     %    of elements. first echo 1-3 of profile_1, angle_1. then echo 1-3 of angle_2. etc...
xyz_check2    = xyz_check( ~isnan(r_backshaped), : );                     % remove all the rows where r is zero (no measurement)

xyz_ori = cat(4,x,y,z);                           % Do the same with original x,y,z: put together in single matrix
xyz_ori = reshape(xyz_ori, [],3);                 % make matrix x,y,z, with each column Txn_profiles
xyz_ori = xyz_ori( ~isnan(reshape(x,[],1)), :);   % select rows without nan

max(abs(xyz_check2 - xyz_ori))  % if different, check if both are in m
    %}

    %% Determine max range, filter
    % plot_interval = 180; % One line for every nth cloud. 180 cloulds=1800 seconds = 1line/30 minutes

    max_r      = max(radius_lidar,[],4);            % Max of the echos, larges distance
%     r_filtered = max(max_r);   % average over clouds
    
    % 99% value
    b = sort(max_r);

    c = floor(sum(~isnan(b)) * 0.99);
    c( ~c) = 1; % if no non-nans, index  is 1 (use first value)
    d = c+[0:719]*size(b,1);
    r_filtered = b(d);
    % filter ----------------------------------------------------------------------------------------------------------------------
    % input filtering
    hoek_use_1 = -90;  % dont plot angles between -180° and -90°
    hoek_use_2 = 105;  % dont plot angles between 105° and 180°
    hoek_frame_skip =   [-17 -15.5; -14 -11; -2.5 -0.5];  % skip the angles -17°to-15.5°; -14°to-11°and -2.5°to-0.5°: frame visible. (for some lines -13.5 and/or -11.5 would be enough)
    x_skip_area = 6;   % delete every point that is below x=6 and (simultaneously) above z=-1m.
    z_skip_area = -1;
    x_plot_max = 12.5; % skip everything past x=12.5m, beyond dune crest (irrespective of z)
    z_plot_max = 0.5; % skip everything above z=0.5m, above dune crest (irrespective of x)

    % filter above lidar
    r_filtered(angle_laserbeam<hoek_use_1)  = nan; % remove around lidar
    r_filtered(angle_laserbeam> hoek_use_2) = nan; % remove around lidar

    % filter frame
    for n=1:3
        r_filtered( (angle_laserbeam>hoek_frame_skip(n,1)) & (angle_laserbeam<hoek_frame_skip(n,2)) ) = nan;
    end

    [x_frame_line, y_frame_line, z_frame_line] = sph2cart(deg2rad(hoek_frame_skip),deg2rad(repmat(angle_laserplane(1,1),3,2)),5*ones(3,2)); % coordinate to plot line of skipped area. in lidar coordinates

    % convert to x, y, z
    [x_lidar,y_lidar,z_lidar] = sph2cart( deg2rad(angle_laserbeam), deg2rad(angle_laserplane), r_filtered ); % convert to xyz coordinates, in mm

    % convert to local coordinates (rough)
    z_local = -x_lidar;
    x_local = y_lidar;
    z_frame_line = - x_frame_line;
    x_frame_line = y_frame_line;

    % filter around lidar (above -1 m, seaward of dune)
    z_local((x_local<x_skip_area) & (z_local>z_skip_area) ) = nan; % remove around lidar

    % filter above crest
    z_local( x_local>x_plot_max) = nan; % remove above crest
    z_local( z_local>z_plot_max ) = nan; % remove above crest

    x_local = x_local'; % column per point cloud
    z_local = z_local';

    %%  Make plot of beach profile every 30 minutes. Use minimum elevation per block to filter out waves   ---------------------------
    plot(x_local,z_local,'.','Color',lijn_kleuren(n_plot,:))  % plot lines

    % % plot patches of skipped areas
%     fill([-50 6 6 -50],[-1 -1 0.5 0.5],[0 0 0],'FaceAlpha',0.1)  % around lidar, delete points
%     fill([12.5 30 30 12.5],[-6 -6 3 3],[0 0 0],'FaceAlpha',0.1)  % landward of crest
%     fill([-50 30 30 -50],[0.5 0.5 3 3],[0 0 0],'FaceAlpha',0.1)  % above crest
%     fill([x_frame_line(1,:),0],[z_frame_line(1,:),0],[0 0 0],'FaceAlpha',0.1)  % frame
%     fill([x_frame_line(2,:),0],[z_frame_line(2,:),0],[0 0 0],'FaceAlpha',0.1)  % frame
%     fill([x_frame_line(3,:),0],[z_frame_line(3,:),0],[0 0 0],'FaceAlpha',0.1)  % frame
%     fill(100*[sind(hoek_use_1) 0 0],100*[cosd(hoek_use_1) 0 1],[1 0 0],'FaceAlpha',0.1)  % sky angle, don't load?
%     fill(100*[sind(hoek_use_2) 0 0],100*[-cosd(hoek_use_2) 0 1],[1 0 0],'FaceAlpha',0.1) % sky angle, don't load?

end
    plot(0,0,'*r') % plot lidar location

% formatting
xlim([-2 12])
ylim([-3.5 0.5])
xlabel('x tov lidar [m]')
ylabel('z tov lidar [m]')
title('Storm 1, Lidar 1')

c = colorbar;
for n = 1 : length(c.TickLabels)   % for every tick label of the colorbar
    n_tekst = n_profile_plot(1) + (n_profile_plot(end)-n_profile_plot(1))/(length(c.TickLabels)-1)*(n-1); % use the nth text string for the date
    c.TickLabels{n} = t_string{round(n_tekst)}(1:16);  % set this text string to the tick label. letter 1:16: 'yyyy-mm-dd hh:mm', i.e. skip seconds
end
ax=gca;
ax.Position(1)=ax.Position(1)-0.03; % move axis to left, to make space for date label
