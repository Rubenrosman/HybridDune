% Try plotting profiles
% NOTES 
%{
- carefull with plotting long durations: you may run out of memory, and need to load data per block instead of for the full duration.
  The current random 6hr period requires 6GB of RAM
- when changing the lidar or storm, change which angles the frame is (line 94), possibly date (line 31-36), 
- Note: this is a random example. In some cases, the 1-minute minimum elevation still shows waves instead of beach

Important lines
- Line 52: Plot file 11, i.e. storm 2, lidar 4
- Line 28: corresponding netcdf file
- Line 33,34 (for storm 2): Date-time to plot
- Line 85: plotting interval: line with min depth per (here) 240 point clouds, so 240/4hz=1min
%}

%% load data
clear, clc
map_files = {'O:\HybridDune experiment\2024-12-18 to 2024-12-20, Storm 1\Lidars\20241220_LiDAR1\storm1_lidar1_polar_10sInterval.nc'
             'O:\HybridDune experiment\2024-12-22, Storm 2\Lidars\20241223_LiDAR1\storm2_lidar1_polar_10sInterval.nc'
             'O:\HybridDune experiment\2025-01-01, Storm 3\Lidars\20250102_LiDAR1\storm3_lidar1_polar_10sInterval.nc'
             'O:\HybridDune experiment\2024-12-18 to 2024-12-20, Storm 1\Lidars\20241220_LiDAR2\storm1_lidar2_polar_10sInterval.nc'
             'O:\HybridDune experiment\2024-12-22, Storm 2\Lidars\20241223_LiDAR2\storm2_lidar2_polar_10sInterval.nc'
             'O:\HybridDune experiment\2025-01-01, Storm 3\Lidars\20250102_LiDAR2\storm3_lidar2_polar_10sInterval.nc'
             'O:\HybridDune experiment\2024-12-18 to 2024-12-20, Storm 1\Lidars\20241220_LiDAR3\storm1_lidar3_polar_10sInterval.nc'
             'O:\HybridDune experiment\2024-12-22, Storm 2\Lidars\20241223_LiDAR3\storm2_lidar3_polar_10sInterval.nc'
             'O:\HybridDune experiment\2025-01-01, Storm 3\Lidars\20250102_LiDAR3\storm3_lidar3_polar_10sInterval.nc'
             'O:\HybridDune experiment\2024-12-18 to 2024-12-20, Storm 1\Lidars\20241220_LiDAR4\storm1_lidar4_polar_10sInterval.nc'
             'O:\HybridDune experiment\2024-12-22, Storm 2\Lidars\20241223_LiDAR4\storm2_lidar4_polar.nc'   % <<------------
             'O:\HybridDune experiment\2025-01-01, Storm 3\Lidars\20250102_LiDAR4\storm3_lidar4_polar_10sInterval.nc'} ;

t1_start = datetime(2024,12,18,23,45,0);
t1_end = datetime(2024,12,20,0,15,0);
t2_start = datetime(2024,12,22,17,00,0);   % storm 2
t2_end = datetime(2024,12,22,23,00,0);     %
t3_start = datetime(2025,1,1,09,45,0);
t3_end = datetime(2025,1,2,10,15,0);
t_storm = [t1_start t2_start t3_start; t1_end t2_end t3_end];
t_storm = datenum(t_storm);
titles_lidar = {'S1 Dike-in-dune'; 'S2 Sandy dune'; 'S3 Dike'; 'S4 Wall-in-dune'};
i_lidar = [1 1 1 2 2 2 3 3 3 4 4 4];
i_storm = [1 2 3 1 2 3 1 2 3 1 2 3];

x_constructions_local = [803.112 791.695 nan     nan     nan
                         nan     nan     nan     nan     nan
                         803.022 791.663 nan     nan     nan 
                         800.033 802.473 802.473 800.033 800.033];
z_constructions_local = [4.78    1.071   nan     nan     nan
                         nan     nan     nan     nan     nan
                         4.564   0.847   nan     nan     nan 
                         1.488   1.488   4.088   4.088   1.488];

for n_file = 11%[2 4 5 11]
    n_lidar = i_lidar(n_file);
    n_storm =i_storm(n_file);
    map_file = map_files{n_file};
    netcdf_info = ncinfo(map_file);
    
    i_max = netcdf_info.Dimensions(1).Length;
    t_num =   ncread(map_file,'time_num',1,i_max,1); % load all times from netcdf

    i1 = find(t_num>t_storm(1,n_storm),1); % point cloud index of start of storm
    i2 = find(t_num>t_storm(2,n_storm),1); % point cloud index of end of storm
%     if isempty(i2) % if measurement already stopped at chosen end of storm
%         i2 = i_max;
%         i1 = i_max - 24*360 - 1;

    i_count = i2-i1+1;
    n_start  = [i1, 1,   3, 1];    % first pointcloud, first point in profile, third profile
    n_count  = [i_count, 720, 1, 3];  % download over how many obs/gaps
    n_stride = [1,   1,   1, 1];  % with which gap width
    
    radius_lidar   = single( ncread(map_file,'radius_lidar',n_start,n_count,n_stride) )/1000; % radius, uint16
    profile_angle  = single( ncread(map_file,'profile_angle', n_start(3),1) )'; % profile angle, 1x16
    ini_beam_angle = single( ncread(map_file,'ini_beam_angle',n_start(3),1) )'; % profile angle, 1x16
    t_string       = ncread(map_file,'time_string',n_start(1),n_count(1),n_stride(1));
    t_num          = ncread(map_file,'time_num'   ,n_start(1),n_count(1),n_stride(1)); % load only times of point clouds used
    
    %%
    angle_laserplane = reshape(profile_angle,1,1,[]);  % in third dimension: third dim om radius is the proifles
    angle_laserbeam  = reshape(ini_beam_angle,1,1,[]); % idem
    angle_laserplane = repmat(angle_laserplane, n_count(1), 720, 1);  % repeat: make T x 720 x n_profiles x 3 
    angle_laserbeam  = repmat(angle_laserbeam, n_count(1), 720, 1) + [0:0.5:359.9];
    
    %% Determine max range, filter
    plot_interval = 240; % One line for every nth cloud. 
    
    max_r      = max(radius_lidar,[],4);            % Max of the echos, largest distance
    r_filtered = movmax(max_r, plot_interval);   % max over clouds
    
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
    
    [x_frame_line, y_frame_line, z_frame_line] = sph2cart(deg2rad(hoek_frame_skip),deg2rad(angle_laserplane(1:3,1:2)),5*ones(3,2)); % coordinate to plot line of skipped area. in lidar coordinates
    
    % convert to x, y, z
    [x_lidar,y_lidar,z_lidar] = sph2cart( deg2rad(angle_laserbeam), deg2rad(angle_laserplane), r_filtered ); % convert to xyz coordinates, in mm
    
    % filter around lidar (above -1 m, seaward of dune)
    x_lidar((y_lidar<x_skip_area) & (-x_lidar > z_skip_area) ) = nan; % remove around lidar
    
    % filter above crest. temp code, in lidar coordinates
    x_lidar( y_lidar>x_plot_max) = nan; % remove above crest
    x_lidar( -x_lidar>z_plot_max ) = nan; % remove above crest
    
    % convert to local coordinates new
    xyz_lidar = [permute(x_lidar,[2 3 1]),  permute(y_lidar,[2 3 1]),  permute(z_lidar,[2 3 1])];
    xyz_local = convert_lidar_coordinates(xyz_lidar,n_lidar, n_storm);
    lidar_position_local = convert_lidar_coordinates([0 0 0],n_lidar, n_storm);
    x_local   = squeeze(xyz_local(:,1,:));
    z_local   = squeeze(xyz_local(:,3,:));
    
    % % convert to local coordinates (rough)
    % z_local = -x_lidar;
    % x_local = y_lidar;
    % z_frame_line = - x_frame_line;
    % x_frame_line = y_frame_line;
    
    % x_local = x_local'; % column per point cloud
    % z_local = z_local';
    
    %%  Make plot of beach profile every 30 minutes. Use minimum elevation per block to filter out waves   ---------------------------
    fig=figure(4); clf(4), hold on
    n_profile_plot = plot_interval : plot_interval : i_count-plot_interval;
    
    % set colors of lines
    colorcet('R1','N',length(n_profile_plot)) % set colors of plot, to ensure the colorbar hasthe same colorscale as the lines
    lijn_kleuren = colorcet('R1','N',length(n_profile_plot)); % define colors of the lines. Set a rainbow color scale with as many colors as no. of lines

    % plot lidar, constructions
    plot(lidar_position_local(1),lidar_position_local(3),'*r') % plot lidar location
    plot(x_constructions_local(n_lidar,:),z_constructions_local(n_lidar,:),'k','LineWidth',1.5)

    for n = 1: length(n_profile_plot)
        n_profile = n_profile_plot(n);
        plot(x_local(:,n_profile),z_local(:,n_profile,:),'-','Color',lijn_kleuren(n,:))  % plot lines
        date_temp = datetime(t_num(n_profile_plot(n)), 'ConvertFrom','datenum');
        legend_string{n} = string(date_temp,'dd MMM HH:mm');
    end
    
    
    % formatting
    if n_lidar ==2
        xlim([788 801])
    else
        xlim([796 801])
    end
    ylim([0.5 5])
    xlabel('x [m]','FontSize',12)
    ylabel('Bed level [m NAP]','FontSize',12)
    title([titles_lidar{n_lidar}, ', Storm ',num2str(n_storm)])
    title('S4,storm 2, individual point clouds, 5 min window')
    % L = legend(legend_string, 'Location','SE')
    % L.Title.String = 'Legend';
    c = colorbar;
    c.Ticks = 0 : 1/8:1;
    for n = 1 : length(c.TickLabels)   % for every tick label of the colorbar
        n_tekst = round( n_profile_plot(1) + (n_profile_plot(end)-n_profile_plot(1))/(length(c.TickLabels)-1)*(n-1) ); % use the nth text string for the date
        date_temp = datetime(t_num(n_tekst), 'ConvertFrom','datenum');
        c.TickLabels{n} = string(date_temp,'dd-MMM HH:mm');
    end
    c.YDir = 'reverse';
    ax=gca;
    ax.Position(1)=ax.Position(1)-0.02; % move axis to left, to make space for date label
    if n_lidar == 2
        ax.Position(3)=0.577; % x bereik kleiner, dus grafiek smaller
    end

    ax.Color='none';
end

%%
% z(x<799.5) = nan;
% z(x>800.5) = nan;
% size(min(z))
% 
% t = datetime(t_num(n_profile_plot), 'ConvertFrom','datenum');
% figure(5)
% plot(t,min(z))
% ylabel('minimum observed elevation between x=799.5 and x=800.5 (just in front of container) [m NAP]')