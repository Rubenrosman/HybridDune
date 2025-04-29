clear, clc

% Note: start has some shenenigans for plotting other storms/lidars (carry-over from other scripts). 
map_file =              'C:\Matlab\hybrid_dunes\storm2_lidar1_polar_10sInterval.nc';
% titles = {'Storm 1, Lidar 1'; 'Storm 2, Lidar 1'; 'Storm 1, Lidar 2'; 'Storm 2, Lidar 2'; 'Storm 1, Lidar 3'; 'Storm 2, Lidar 3'; 'Storm 1, Lidar 4'; 'Storm 2, Lidar 4';};
t1_start = datetime(2024,12,18,23,45,0);
t1_end =   datetime(2024,12,20,0,15,0);
t2_start = datetime(2024,12,20,14,45,0);  % effectively from 15:00:frist value is avg 14:45-15:15
t2_end =   datetime(2024,12,23,15,15,0);  % 3 days later
t3_start = datetime(2025,1,1,09,45,0);
t3_end = datetime(2025,1,2,10,15,0);
t_storm = [t1_start t2_start t3_start; t1_end t2_end t3_end];
t_storm = datenum(t_storm);
titles_lidar = {'S1 Dike-in-dune'; 'S2 Sandy dune'; 'S3 Dike'; 'S4 Wall-in-dune'};
i_lidar = [1 1 1 2 2 2 3 3 3 4 4 4];
i_storm = [1 2 3 1 2 3 1 2 3 1 2 3];

% pre-define coordinates from container, revetment
x_constructions_local = [803.112 791.695 nan     nan     nan
                         nan     nan     nan     nan     nan
                         803.022 791.663 nan     nan     nan 
                         800.033 802.473 802.473 800.033 800.033];
z_constructions_local = [4.70    1.071   nan     nan     nan
                         nan     nan     nan     nan     nan
                         4.564   0.847   nan     nan     nan 
                         1.488   1.488   4.088   4.088   1.488];

for n_file = [2]
    n_lidar = i_lidar(n_file);
    n_storm =i_storm(n_file);

    %% load data ------------------------------------------
    netcdf_info = ncinfo(map_file);
    
    i_max = netcdf_info.Dimensions(1).Length;
    t_num =   ncread(map_file,'time_num',1,i_max,1); % load all times from netcdf

    i1 = find(t_num>t_storm(1,n_storm),1); % point cloud index of start of storm
    i2 = find(t_num>t_storm(2,n_storm),1); % point cloud index of end of storm
%     if isempty(i2) % if measurement already stopped at chosen end of storm
%         i2 = i_max;
%         i1 = i_max - 24*360 - 1;
%     i1=1;
%     i2 = i_max;

    i_count = i2-i1+1;
    n_start  = [i1, 1,   3, 1];    % first pointcloud, first point in profile, third profile
    n_count  = [i_count, 720, 1, 3];  % download over how many obs/gaps
    n_stride = [1,   1,   1, 1];  % with which gap width
    
    % load datafile (cross-shore profile (no3), full file from start to end, selecting of timestep comes later)
    radius_lidar   = single( ncread(map_file,'radius_lidar',n_start,n_count,n_stride) )/1000; % radius, uint16. Size T x 720 angles x 1 profile x 3 echos
    profile_angle  = single( ncread(map_file,'profile_angle', n_start(3),1) )'; % profile angle, 1x16
    ini_beam_angle = single( ncread(map_file,'ini_beam_angle',n_start(3),1) )'; % profile angle, 1x16
    t_string       = ncread(map_file,'time_string',n_start(1),n_count(1),n_stride(1));
    t_num          = ncread(map_file,'time_num'   ,n_start(1),n_count(1),n_stride(1)); % load only times of point clouds used
    
    % calculate angles corresponding to the cells of the radius matrix
    angle_laserplane = reshape(profile_angle,1,1,[]);  % in third dimension: third dim om radius is the proifles
    angle_laserbeam  = reshape(ini_beam_angle,1,1,[]); % idem
    angle_laserplane = repmat(angle_laserplane, n_count(1), 720, 1);  % repeat: make T x 720 x n_profiles x 3 
    angle_laserbeam  = repmat(angle_laserbeam, n_count(1), 720, 1) + [0:0.5:359.9];
    
    %% Determine max range, filter
    plot_interval = 180; % One line for every nth cloud. 180 cloulds=1800 seconds = 1line/30 minutes

    max_r      = max(radius_lidar,[],4);            % Max of the echos, largest distance
    r_filtered = movmax(max_r, plot_interval);   % average over clouds
    
    % filter ----------------------------------------------------------------------------------------------------------------------
    % input filtering
    hoek_use_1 = -90;  % dont plot angles between -180° and -90°
    hoek_use_2 = 105;  % dont plot angles between 105° and 180°
              
    % 3 tubes that block view of lidar. for tube 1/2/3, skip respectively between value 1 and 2 / between 3&4 / between 5&6
    hoek_frame_skip =   [-16.9 -13.9; -13.3 -9.9; -2.25 0.5];  % S1, storm 2

    x_skip_area = 6;   % delete every point that is below x=6 and (simultaneously) above z=-1m. 
    z_skip_area = -1; 
    x_plot_max = 12.5; % skip everything past x=12.5m, beyond dune crest (irrespective of z)
    z_plot_max = 0.5; % skip everything above z=0.5m, above dune crest (irrespective of x)
    
    % filter above lidar
    r_filtered(angle_laserbeam<hoek_use_1)  = nan; % remove around lidar
    r_filtered(angle_laserbeam> hoek_use_2) = nan; % remove around lidar
    
    % filter frame
    for n=1:3 % filter 3 tubes of frame
        r_filtered( (angle_laserbeam>hoek_frame_skip(n,1)) & (angle_laserbeam<hoek_frame_skip(n,2)) ) = nan;
    end
    r_filtered( (angle_laserbeam>-50) & (angle_laserbeam<-40 ) )= nan;

    [x_frame_line, y_frame_line, z_frame_line] = sph2cart(deg2rad(hoek_frame_skip),deg2rad(angle_laserplane(1:3,1:2)),5*ones(3,2)); % coordinate to plot line of skipped area. in lidar coordinates
    
    % convert to x, y, z
    [x_lidar,y_lidar,z_lidar] = sph2cart( deg2rad(angle_laserbeam), deg2rad(angle_laserplane), r_filtered ); % convert to xyz coordinates, in mm
    
    % filter around lidar (above -1 m, seaward of dune)
    x_lidar((y_lidar<x_skip_area) & (-x_lidar > z_skip_area) ) = nan; % remove around lidar

    % filter above crest. temp code, in lidar coordinates
    x_lidar( y_lidar>x_plot_max) = nan; % remove above crest
    x_lidar( -x_lidar>z_plot_max ) = nan; % remove above crest
    
    % convert to local coordinates new
    % convert with function, that wants xyz-matrix in size (Tx3 or Tx3xN).
    xyz_lidar = [permute(x_lidar,[2 3 1]),  permute(y_lidar,[2 3 1]),  permute(z_lidar,[2 3 1])]; % reshape x_lidar from T x 720angles to 720x1xT, and stack with y and z)
    xyz_local = convert_lidar_coordinates(xyz_lidar,n_lidar, n_storm); % function to convert to local coordinates (corner container as 200,800)
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
    fig=figure(2); clf(2), hold on
    n_profile_plot = 91 : plot_interval : i_count; % lines as 30min avg. Start after 15minutes, so the 15minx60s/10sinterval=90th line
    
    % set colors of lines
    colorcet('R1','N',length(n_profile_plot)) % set colors of plot, to ensure the colorbar hasthe same colorscale as the lines
    lijn_kleuren = colorcet('R1','N',length(n_profile_plot)); % define colors of the lines. Set a rainbow color scale with as many colors as no. of lines

    % plot lidar, constructions
%     plot(lidar_position_local(1),lidar_position_local(3),'*r') % plot lidar location
    plot(x_constructions_local(n_lidar,:),z_constructions_local(n_lidar,:),'k','LineWidth',2)

    for n = 1: length(n_profile_plot)
        n_profile = n_profile_plot(n);
        xy_plot = [x_local(:,n_profile),z_local(:,n_profile,:)];
        xy_plot = removegap_maxdif(xy_plot,[0.5 0.2]); % remove gap (NaN) if dx_gap<=0.5m AND dz_gap<=0.2m
        plot(xy_plot(:,1),xy_plot(:,2),'-','Color',lijn_kleuren(n,:),'LineWidth',1.5)  % plot lines
        date_temp = datetime(t_num(n_profile_plot(n)), 'ConvertFrom','datenum');
        legend_string{n} = string(date_temp,'dd MMM HH:mm');
    end
    
    % formatting
    axis equal
    xlim([786.2 804])
    ylim([0.8 5])
    xlabel('x [m]','FontSize',12)
    ylab = ylabel('Elevation [m NAP]','FontSize',12);
    
    % L = legend(legend_string, 'Location','SE')
    % L.Title.String = 'Legend';

     % add and format colorbar
    c = colorbar;
    c.Ticks = 0 : 1/6:1; % 7 ticks/labels (3days, so label every 12 hours)
    for n = 1 : length(c.TickLabels)   % for every tick label of the colorbar
        n_tekst = n_profile_plot(1) + (n_profile_plot(end)-n_profile_plot(1))/(length(c.TickLabels)-1)*(n-1); % use the nth text string for the date
        date_temp = datetime(t_num(n_tekst), 'ConvertFrom','datenum');
        c.TickLabels{n} = string(date_temp,'dd-MMM HH:mm');
    end
    c.YDir = 'reverse';
    ax=gca;
    ax.Position(1)=ax.Position(1)-0.02; % move axis to left, to make space for date label
    ax.Position(2)=ax.Position(2)+0.01; % move axis up

    grid on
    ax.Color='none';
    ax.GridAlpha = 0.5;
%     ax.YRuler.TickLabelGapOffset = -22;
%     ylab.Position(1) = ylab.Position(1)-0.5
fig.Position=[100 100 950 360];

end

%% export_fig(fig,'C:\Users\dpoppema\OneDrive - Delft University of Technology\PostDoc Hybrid dunes\Conferences\2025 - NCK\poster\plot4.png','-m6')
locatie = ['C:\Users\dpoppema\OneDrive - Delft University of Technology\PostDoc Hybrid dunes\Conferences\2025 - Coastal Dynamics\presentatie\figures\'];
naam    = ['S1 Storm2 30min bedlevel dynamics.png'];

% exportgraphics(fig,[locatie,naam],'Resolution','300')

function matrix_filled = removegap_maxdif(matrix_Nx2,max_gap_1x2) 
% remove rows with nans from matrix for gaps up to size max_gap, i.e. if the values before and after the gap differ at most by max_gap_1x2 
%{
input: 
- a Nx2 xy matrix, with rows of NaN's (i.e. if x=nan, then always y=nan). 
- the max gap size, same unit as Nx2

example: 
xy = [1   1
      nan nan
      nan nan
      2   1.5], in meters 
The x-gap is 1m, the y-gap is 0.5 meters. So removegap_maxdif removes the nan rows if the first value of max_gap_1x2 at least 1 AND the second value at least 0.5
%}
    M = isnan(matrix_Nx2);
    
    % determine consequetive rows with nan
    gap = zeros(length(M),size(matrix_Nx2,2)); % in unit of matrix_Nx2

    for i = 2 : length(M)  % count consequetive nan
        if M(i) % if nan
            i_start = find(~M(1:i-1,1),1,'last'); % index of last value that wasn't nan
            i_end = find(~M(i+1:end,1),1)+i;      % index of first next value that isn't nan
            if isempty(i_start) || isempty(i_end) % if no previous/next value, e.g. matrix starts with nan
                gap(i,:)=nan;                     % then gap size is nan
            else
                gap(i,:) = abs( matrix_Nx2(i_end,:) - matrix_Nx2(i_start,:) ); %calculate difference between previous and next value, i.e. 2 meter
            end
        end
    end
    
    % fill gap
    keep_rows = any( ~M | (gap >max_gap_1x2) ,2);
    matrix_filled = matrix_Nx2(keep_rows,:);
end