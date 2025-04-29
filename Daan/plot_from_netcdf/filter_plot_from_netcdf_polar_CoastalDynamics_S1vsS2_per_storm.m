% Plot profiles Coastal Dynamics presenation

%{
1 define file names
2 define dates, times
3 pre-define plot
4 load data
5 select desired date, time
6 filter and georeference
7 plot
%}
clear, clc
for n_plot = 8  % make plot for every storm. n=1-4: only S2. n=5-8: also S1. 1: initial. 2/3/4: after storm1/2/3
    n_geplot=0;
    %% 1) Define file names
    clc
    % map-files: rows: n_lidar. cols: n_storm
    map_files = {'O:\HybridDune experiment\storm1_lidar1_polar_10sInterval.nc',...
                 'O:\HybridDune experiment\storm2_lidar1_polar_10sInterval.nc',...
                 'O:\HybridDune experiment\storm3_lidar1_polar_10sInterval.nc';
                 'O:\HybridDune experiment\storm1_lidar2_polar_10sInterval.nc',... % lidar 2
                 'O:\HybridDune experiment\storm2_lidar2_polar_10sInterval.nc',...
                 'O:\HybridDune experiment\storm3_lidar2_polar_10sInterval.nc';
                 'O:\HybridDune experiment\storm1_lidar3_polar_10sInterval.nc',... % lidar 3
                 'O:\HybridDune experiment\storm2_lidar3_polar_10sInterval.nc',...
                 'O:\HybridDune experiment\storm3_lidar3_polar_10sInterval.nc';
                 'O:\HybridDune experiment\storm1_lidar4_polar_10sInterval.nc',... % lidar 4
                 'O:\HybridDune experiment\storm2_lidar4_polar_10sInterval.nc',...
                 'O:\HybridDune experiment\storm3_lidar4_polar_10sInterval.nc'} ;
    
    %% 2) Define dates, times
    datetime_high = {[2024,12,19, 5, 0,0; % storm 1
                      2024,12,19,17,30,0]; %              peak
                     [2024,12,22, 7,30,0; % storm 2      peak
                      2024,12,22,20, 0,0];
                     [2025, 1, 1,15, 0,0; % storm 3      peak
                      2025, 1, 2, 4, 0,0;
                      2025, 1, 2,16,30,0]};
    
    datetime_low = {[ 2024,12,18,23, 0,0; % storm 1      
                      2024,12,19,11,15,0; %              before peak
                      2024,12,19,23,30,0];
                     [2024,12,22, 1,30,0; % storm 2      before peak
                      2024,12,22, 13,45,0;
                      2024,12,23, 2, 0,0];
                     [2025, 1, 1, 9, 0,0; % storm 3      before peak
                      2025, 1, 1,21,15,0;
                      2025, 1, 2, 10,15,0;
    %                  2025, 1, 2,22,30,0;... % skip, not in lidar
                      ]};
    
    x_constructions_local = [803.112 791.695 nan     nan     nan
                             nan     nan     nan     nan     nan
                             803.022 791.663 nan     nan     nan 
                             800.033 802.473 802.473 800.033 800.033];
    z_constructions_local = [4.70    1.071   nan     nan     nan
                             nan     nan     nan     nan     nan
                             4.564   0.847   nan     nan     nan 
                             1.488   1.488   4.088   4.088   1.488];
    %% 3) pre-define plot
    title_string_all = {'S1 Dike-in-dune';'S2 Sandy dune';'S3 Dike';'S4 Wall-in-dune'};
    fig=figure(5); clf(5), hold on
    fig.Position=[100 100 950 360];
    % define colors of the lines. Set a rainbow color scale with as many colors as no. of lines
    % colorcet('R1','N',9) % set colors of plot, to ensure the colorbar hasthe same colorscale as the lines
    % lijn_kleuren = colorcet('R1','N',19);          % version for 4 storms, 3 tides each
    % n_all = [1 2 3   6 7 8   11 12 13   16 17 18]
    lijn_kleuren = zeros(9,3);
    temp = colorcet('R1','N',4);          % version for 3 storms, 3 tides each
    lijn_kleuren(1,:) = temp(1,:);
    lijn_kleuren(3,:) = temp(2,:);
    lijn_kleuren(6,:) = temp(3,:);
    lijn_kleuren(9,:) = temp(4,:);
    % plot(x_constructions_local(n_lidar,:),z_constructions_local(n_lidar,:),'k','LineWidth',2)
    
    for n_lidar = [2 1]%1:2%:4
    
    %% 4) load data
        for n_storm = 1:3
            map_file = map_files{n_lidar,n_storm};  % select netcdf file
            netcdf_info = ncinfo(map_file);         % read number of obs from netcdf
            i_max = netcdf_info.Dimensions(1).Length;
            
            % select which parts of file to load (profile 3, all times and angles
            n_start  = [1, 1,   3, 1];    % first pointcloud, first point in profile, third profile
            n_count  = [i_max, 720, 1, 3];  % download over how many obs/gaps. 
            n_stride = [1,   1,   1, 1];  % with which gap width
            
            % load datafile (cross-shore profile (no3), full file from start to end, selecting of timestep comes later)
            radius_lidar   = single( ncread(map_file,'radius_lidar',n_start,n_count,n_stride) )/1000; % radius, uint16. Size T x 720 angles x 1 profile x 3 echos
            profile_angle  = single( ncread(map_file,'profile_angle', n_start(3),1) )'; % profile angle, 1x16
            ini_beam_angle = single( ncread(map_file,'ini_beam_angle',n_start(3),1) )'; % profile angle, 1x16
            t_string       = ncread(map_file,'time_string',n_start(1),n_count(1),n_stride(1));
            t_num          = ncread(map_file,'time_num'   ,n_start(1),n_count(1),n_stride(1)); % load only times of point clouds used
    
    %% 5) Select desired date, time
            datetime_low_num = datetime(datetime_low{n_storm}); % make datetime array: low tides of storm n
            datetime_low_num = datenum(datetime_low_num ); % convert to num array
            clear i_low_tide
            for n = 1:length(datetime_low_num)
                i_low_tide(n) = find(t_num>datetime_low_num(n),1); % point cloud index of low tide n of storm
            end
            
            % select data from these times
            radius_lidar = radius_lidar(i_low_tide,:,:,:); % select timestep. size T (i.e. few low tides) x 720 x 1 profile x 3 echos
            radius_lidar = max(radius_lidar,[],4); % % Max of the echos, largest distance, to skip water echos. Resulting size: low tides (rows) x 720 cols
    
    %% 6 Filter and georeference
            % define angles at which the lidar sees the frame. These will be skipped in plotting.
            hoek_frame_skip(:,:,1,1) =   [-17.4 -14.9; -13.8 -10.9; -2.75 0];  % S1, storm 1
            hoek_frame_skip(:,:,1,2) =   [-16.9 -13.9; -13.3 -9.9; -2.25 0.5];  % S1, storm 2
            hoek_frame_skip(:,:,1,3) =   [-17.4 -14.4; -13.8 -10.9; -2.75 0];  % S1, storm 3
            
            hoek_frame_skip(:,:,2,1) =   [-21.2 -18.2; -14.3 -11.3; -2.3 0.7];  % S2, storm 1
            hoek_frame_skip(:,:,2,2) =   [-20.2 -17.7; -13.8 -10.8; -1.3 1.7];  % S2, storm 2
            hoek_frame_skip(:,:,2,3) =   [-21.2 -18.2; -14.3 -11.8; -2.3 0.7];  % S2, storm 3
            
            hoek_frame_skip(:,:,3,1) =   [-17.0 -14.5; -12.5 -9.5; 1.0 4.0];  % S3, storm 1
            hoek_frame_skip(:,:,3,2) =   [-17.0 -14.5; -12.5 -9.5; 1.0 4.0];  % S3, storm 2
            hoek_frame_skip(:,:,3,3) =   [-17.5 -14.5; -12.5 -10.0; 0.5 4.0];  % S3, storm 3
            
            hoek_frame_skip(:,:,4,1) =   [-19.25 -16.25; -11.25 -8.75; -2.75 0.25];  % S4, storm 1
            hoek_frame_skip(:,:,4,2) =   [-19.25 -16.25; -11.25 -8.75; -2.75 0.25];  % S4, storm 2
            hoek_frame_skip(:,:,4,3) =   [-19.75 -16.75; -11.75 -8.75; -3.25 -0.25];  % S4, storm 3
            hoek_frame_skip = hoek_frame_skip(:,:,n_lidar,n_storm);
    
            % calculate angles corresponding to the cells of the radius matrix
            n_timesteps = length(i_low_tide); % number of timesteps to plot
            angle_laserplane = reshape(profile_angle,1,1,[]);  % in third dimension: third dim of radius is the profiles (only relevant if multiple profiles used, otherwise single value). 
            angle_laserbeam  = reshape(ini_beam_angle,1,1,[]); % idem
            angle_laserplane = repmat(angle_laserplane, n_timesteps, 720, 1);  % repeat: make T x 720 
            angle_laserbeam  = repmat(angle_laserbeam, n_timesteps, 720, 1) + [0:0.5:359.9];
    
            % filter in polar coordinates  -----------------------------------------------------------------------------------------------------
            % input filtering. NOTE: FILTERING HERE TO FIND BEACH PROFILES. SOME WAVES FILTERED OUT
            hoek_use_1 = -90;  % dont plot angles between -180° and -90° (between vertical and horizontal seaward, assumed to be rain/noise)
            hoek_use_2 = 105;  % dont plot angles between 105° and 180°  (between vertical and slightly upward landward, idem)
            x_skip_area = 6;   % delete every point that is below x=6 (6 m landward of lidar) and (simultaneously, next line) ...
            z_skip_area = -1;    % % above z=-1m (1 m below the lidar or higher)
            x_plot_max = 12.5; % skip everything past x=12.5m, beyond dune crest (irrespective of z)
            z_plot_max = 0.5; % skip everything above z=0.5m, above dune crest (irrespective of x)
    
            % filter above lidar
            radius_lidar(angle_laserbeam<hoek_use_1)  = nan; 
            radius_lidar(angle_laserbeam> hoek_use_2) = nan; 
    
            % filter frame
            for n=1:3
                radius_lidar( (angle_laserbeam>hoek_frame_skip(n,1)) & (angle_laserbeam<hoek_frame_skip(n,2)) ) = nan;
            end
    
            % convert polar to x, y, z ------------
            % NB: lidar coordinates, so lidar still as origin (0,0,0)
            [x_lidar,y_lidar,z_lidar] = sph2cart( deg2rad(angle_laserbeam), deg2rad(angle_laserplane), radius_lidar ); % convert to xyz coordinates, in mm
            [x_frame_line, y_frame_line, z_frame_line] = sph2cart(deg2rad(hoek_frame_skip),deg2rad(angle_laserplane(1:3,1:2)),5*ones(3,2)); % coordinate to plot line of skipped area. in lidar coordinates
    
            % filter in cartesian lidar coordinates------------
            % filter around lidar (above -1 m, seaward of dune)
            x_lidar((y_lidar<x_skip_area) & (-x_lidar > z_skip_area) ) = nan; % remove around lidar
    
            % convert lidar coordinates to local coordinates
            % convert with function, that wants xyz-matrix in size (Tx3 or Tx3xN).
            xyz_lidar = [permute(x_lidar,[2 3 1]),  permute(y_lidar,[2 3 1]),  permute(z_lidar,[2 3 1])]; % reshape x_lidar from 3tides x 720angles to 720x1x3, and stack with y and z)
            xyz_local = convert_lidar_coordinates(xyz_lidar,n_lidar, n_storm); % function to convert to local coordinates (corner container as 200,800)
            x_local   = squeeze(xyz_local(:,1,:)); 
            z_local   = squeeze(xyz_local(:,3,:));
    
            % temporarily filter straight line up above container (paaltje?)
            if n_storm == 2 && n_lidar == 4
                z_local((z_local>4.76) & [1 1 0] ) = nan;
                warning('manual filtering S4, storm2, tide1&2')
            end
%             % convert to local coordinates (rough)
%             z_local = -x_lidar;
%             x_local = y_lidar;
%             z_frame_line = - x_frame_line;
%             x_frame_line = y_frame_line;
%             
%             x_local = x_local'; % column per point cloud
%             z_local = z_local';

    %% 7 plot) 
            for i_profile = 1:3 % code (from older script) to plot 3 profiles per storm. Initial, first low tide, second low tide
                if (n_storm == 1 && i_profile==1) | i_profile == 3 % this time only plot 4 lines. The initial dune profile, and final profile (second low tide) per storm
                    n_geplot=n_geplot+1;
    
                    n_kleur = (n_storm-1)*3 + i_profile; % recall color, set as start of script
                    xy_plot = [x_local(:,i_profile),z_local(:,i_profile,:)];
                    xy_plot = removegap_maxdif(xy_plot,[0.5 0.2]); % remove gap (NaN) if dx_gap<=0.5m AND dz_gap<=0.2m
                    flag_visible = n_plot>=n_geplot; % to make separate plots of different stages
                    line_width = 0.75 + 0.75*(n_plot == n_geplot) + 0.5*(n_plot -4 == n_geplot); % different line width per plot
                    if n_lidar == 1
                        p{n_geplot}=plot(xy_plot(:,1),xy_plot(:,2),'--','Color',lijn_kleuren(n_kleur,:),'LineWidth',line_width,'Visible',flag_visible);  % plot lines
                    else
                        p{n_geplot}=plot(xy_plot(:,1),xy_plot(:,2),'-','Color',lijn_kleuren(n_kleur,:),'LineWidth',line_width,'Visible',flag_visible);  % plot lines
                    end
                end
            end
        end
    %     fill([x_frame_line(1,:),0],[z_frame_line(1,:),0],[0 0 0],'FaceAlpha',0.1)  % plot frame. only works when plotting in rough local coordinates, centered around 0
    %     fill([x_frame_line(2,:),0],[z_frame_line(2,:),0],[0 0 0],'FaceAlpha',0.1)  % frame
    %     fill([x_frame_line(3,:),0],[z_frame_line(3,:),0],[0 0 0],'FaceAlpha',0.1)  % frame
    if n_plot<5
        p{9}=plot(x_constructions_local(n_lidar,:),z_constructions_local(n_lidar,:),'k','LineWidth',2,'Visible',0);
    else
        p{9}=plot(x_constructions_local(n_lidar,:),z_constructions_local(n_lidar,:),'k','LineWidth',2,'Visible',1);
    end

    end % ends for n_lidar=... loop
    
        xlabel('Cross-shore distance [m]')
        ylabel('Evelevation [m NAP]')
    %     title('S1 Dike-in-dune vs S2 Sandy dune')
    
        grid on
        
    %     plot(0,0,'*r')
    %     xlim([-1.5 0.5])
    
        axis equal
        xlim([786 804])
        ylim([0 6])   
        legend_text = {'S1, Before storm 1','S1, After storm 1','S1, After storm 2','S1, After storm 3','S1, Revetment',...
                       'S2, Before storm 1','S2, After storm 1','S2, After storm 2','S2, After storm 3'};
    
        L = legend([p{5} p{6} p{7} p{8} p{9} p{1} p{2} p{3} p{4}],legend_text,'Location','NW','NumColumns',2);
        ax=gca;
        L.Position(1)=0.1299; 
        L.Position(2)=0.6222; 
    
    %% save figure
    locatie = ['C:\Users\dpoppema\OneDrive - Delft University of Technology\PostDoc Hybrid dunes\Conferences\2025 - Coastal Dynamics\presentatie\figures\S1 vs S2\'];
    naam    = ['S1vsS2 V',num2str(n_plot),'.png'];
    
%     exportgraphics(fig,[locatie,naam],'Resolution','300') % uncomment to save figure
end % ends for n_plot loop

%%
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