%% Load all raw data: OSSI
clear, clc
map = 'O:\HybridDune experiment\data RBR, OSSI\Ossi data\raw_netcdf\';
t_cal_factor = [1.000044874572203   1.000042922146781   1.000034643947682   1.000025897320181   1.000032925848062   1.000041360211933   1.000036518244910   1.000025038284660];
offset_OSSI = [ 84667   99069   87006   96112   93065   83751   94704   93712]; % based on the second bucket test

for n=1:8
    map_naam = [map, 'OSSI ', num2str(n), '.nc'];
%     a = ncinfo(map_naam);
    p = ncread(map_naam, 'p');                    % load time + attributes of time
    t = ncread(map_naam, 't');                    % load time + attributes of time
    t_atts = ncinfo(map_naam,'t').Attributes;      % Load attributes of time, select unit used
    t = read_time_from_xarray_netcdf(t,t_atts);         % convert to datetime format

    p_all_OSSI{n} = single(p) / 1e4; % pressure in m
    t_all_OSSI{n} = t;
end

%% Load all raw data: RBR
map = 'O:\HybridDune experiment\data RBR, OSSI\copy RBR Udrive series1\raw NetCDF\';
names = {'refP1 RBR4', 'S1P3 RBR5', 'S2P3 RBR1', 'S3P3 RBR6', 'S4P3 RBR2','S1P2 RBR3'};
% offset_RBR        = [-597,  -225,   315,   30,    353,    124];  % old values, calibration period between 12dec 20:00 and 13dec 6:00
offset_RBR        = [-543,  -147,   268,   38     320      64];    % new values, based on 23dec between 18h and 19h, just before the final bucket test
warning('check if RBR offset is for the correct period')
for n=1:6
    map_naam = [map, names{n}, ' raw data - period 1.nc'];
%     a = ncinfo(map_naam);
    p = ncread(map_naam, 'p');                    % load time + attributes of time
    t = ncread(map_naam, 't');                    % load time + attributes of time
    t_atts = ncinfo(map_naam,'t').Attributes;      % Load attributes of time, select unit used
    t = read_time_from_xarray_netcdf(t,t_atts);         % convert to datetime format

    p_all_RBR{n} = ( single(p) + offset_RBR(n) ) / 1e4; % pressure in m
    t_all_RBR{n} = t;
end

%%
clc
% offset_OSSI = [ 84717   99188   87040   96189   93167   83849   94742   93758]; % based on second bucket test
% offset_OSSI = [ 84890   99142   87109   96668   93117   83934   94020   93639]; % based on first bucke test
% offset_OSSI = [ 84667   99069   87006   96112   93065   83751   94704   93712]; % just before the last bucket test

% Set 8 different colors for plotting OSSIs
colors_OSSI = zeros(8,3);
for n=1:6
    colors_OSSI(n,:)=default_color(n);
end
colors_OSSI(7:8,:) = [0 1 1; 1 0 1];

for n=1:8
    OSSI_labels{n} = ['OSSI', num2str(n)];
end
%% Plot all OSSI sensors
figure(1), clf(1), hold on
for n=1:8
    t = t_all_OSSI{n}; % + dt_OSSI(n)/24/3600;
    p = p_all_OSSI{n} *1e4 + offset_OSSI(n);
    i = 1 : 10 : round( length(t) ); % plot every 0.5s 

    plot( t(i), p(i), 'color',colors_OSSI(n,:))  % Plot calibrated data
end
% xlim([t1 t2])
% ylim([10 11])
legend(OSSI_labels)
grid on
ylabel('Pressure [Pa]')

%% Figure of timing at start of test: timing OSSIs correct
t0_end = [datetime(2024,12,16,  8,00,00), datetime(2024,12,16,  9,00,00)   % sensor start
          datetime(2024,12,16,  8,11,00), datetime(2024,12,16,  8,12,00)];   % zoom in on peak during sensor start
titles = {'OSSI sensors at start','zoom of OSSI&RBR sensors at start'};

figure(2), clf(2)
for n_subplot = 1:2
    subplot(2,1,n_subplot), hold on
    for n=1:8 % for every subplot, plot 8 lines (8 OSSI sensors)
        [t, p] = select_t_p(t_all_OSSI,p_all_OSSI, n, t0_end(n_subplot,:));  % select t and p: sensor number and time period
        p = p + + offset_OSSI(n)/1e4;
        plot(t,p,'color',colors_OSSI(n,:))         % Plot calibrated data
    end

    if n_subplot==1
        legend(OSSI_labels)
    else % only plot rbr for second subplot
    for n=1:6 % for every subplot, plot 6 RBR sensors
        [t, p] = select_t_p(t_all_RBR,p_all_RBR, n, t0_end(n_subplot,:));  % select t and p: sensor number and time period
        plot(t,p,'k')
    end
        % mark period
        plot([datetime(2024,12,16, 8,11,22) datetime(2024,12,16, 8,11,22)], [10.2 10.6],'k')
        plot([datetime(2024,12,16, 8,11,30) datetime(2024,12,16, 8,11,30)], [10.2 10.6],'k')
        plot([datetime(2024,12,16, 8,11,22) datetime(2024,12,16, 8,11,30)], [10.2 10.6; 10.2 10.6], 'k')
    end
    
    title(titles{n_subplot})
    grid on
end

%% Calibrate timing
% Select period of RBR data with clear peak, and determine date_time of peak -----------------------------------------
t0_end = [datetime(2024,12,23, 17,52,22), datetime(2024,12,23, 17,52,30) ];  % period of last peak
for n=1:6
    [t, p] = select_t_p(t_all_RBR,p_all_RBR, n, t0_end);  % select t and p: sensor number and time period
    [p_max,i_max] = max(p);
    t_peak(n) = t(i_max);
end
t_peak_RBR = mean(t_peak);

% Select period of OSSI data with the same peak, and determine date_time -----------------------------------------
% NB: a period with 3 peaks in all sensors (OSSI and RBR) is used. For the RBRs, a short period with only the last peak is selected. For the 
% OSSI's a longer period is used,but the last peak is always the highest.
t0_end = [datetime(2024,12,23, 17,51,45), datetime(2024,12,23, 17,52,30)] ;% same period as above, but slightly longer window
for n=1:8
    [t, p] = select_t_p(t_all_OSSI,p_all_OSSI, n, t0_end);  % select t and p: sensor number and time period
    [p_max,i_max] = max(p);
    t_peak_OSSI(n) = t(i_max);
end
dt_OSSI = t_peak_RBR - t_peak_OSSI;

% Plot results -----------------------------------------
t0_end = [datetime(2024,12,23, 17,51,30), datetime(2024,12,23, 17,52,30)   % 
          datetime(2024,12,23, 17,51,30), datetime(2024,12,23, 17,52,30)
          datetime(2024,12,23, 17,52,10), datetime(2024,12,23, 17,52,30)];
titles = {'RBRs', 'OSSIs', 'all sensors, timing OSSIs corrected'};

figure(3), clf(3)
% RBRs
subplot(2,2,1), hold on
for n=1:6 % for every subplot, plot 6 RBR sensors
    [t, p] = select_t_p(t_all_RBR,p_all_RBR, n, t0_end(1,:));  % select t and p: sensor number and time period
    plot(t,p)
end
grid on, legend, title(titles{1})

% OSSIs
subplot(2,2,2), hold on
for n=1:8 % for every subplot, plot 6 RBR sensors
    [t, p] = select_t_p(t_all_OSSI,p_all_OSSI, n, t0_end(2,:));  % select t and p: sensor number and time period
    plot(t, p+ offset_OSSI(n)/1e4)
end
grid on, legend, title(titles{2})

% Subplot 3: RBRs + OSSIs
subplot(2,1,2), hold on
for n=1:6 % for every subplot, plot 6 RBR sensors
    [t, p] = select_t_p(t_all_RBR,p_all_RBR, n, t0_end(3,:));  % select t and p: sensor number and time period
    plot(t,p)
end
for n=1:8 % for every subplot, plot 6 RBR sensors
    [t, p] = select_t_p(t_all_OSSI,p_all_OSSI, n, t0_end(3,:)-dt_OSSI(n));  % select t and p: sensor number and time period
    plot(t+dt_OSSI(n), p+ offset_OSSI(n)/1e4)
end
grid on, legend, title(titles{3})


%% Reload OSSI data, with time shift
map = 'O:\HybridDune experiment\data RBR, OSSI\Ossi data\raw_netcdf\';

for n=1:8
    map_naam = [map, 'OSSI ', num2str(n), '.nc'];
%     a = ncinfo(map_naam);
    p = ncread(map_naam, 'p');                    % load time + attributes of time
    t = ncread(map_naam, 't');                    % load time + attributes of time
    t_atts = ncinfo(map_naam,'t').Attributes;     % Load attributes of time, select unit used
    t = read_time_from_xarray_netcdf(t, t_atts, t_cal_factor(n));  % correct OSSIs for time shift

    p_all_OSSI{n} = single(p) / 1e4; % pressure in m
    t_all_OSSI{n} = t;
end

%% Calibrate sensors, plot bucket tests 
% cal_period_t1 = datetime(2024,12,16,9,43,30);  % 'platteau' during first bucket test
% cal_period_t2 = datetime(2024,12,16,9,44,0);
cal_period_t1 = datetime(2024,12,23,19,11,0);  % 'platteau' during second bucket test
cal_period_t2 = datetime(2024,12,23,19,25,0);
% cal_period_t1 = datetime(2024,12,23,18,0,0);      % just before bucket test 2
% cal_period_t2 = datetime(2024,12,23,19,0,0);
% cal_period_t1 = datetime(2024,12,16,20,0,0);      % low tide 1
% cal_period_t2 = datetime(2024,12,17,0,0,0);

% Calibrate OSSIs ---------------------------------------
% Load data RBR, determine avg pressure
for n=1:6
    t = t_all_RBR{n};
    p = p_all_RBR{n};

    % determine p during calibration period
    i1 = find(t>cal_period_t1,1);
    i2 = find(t>cal_period_t2,1);
    p_mean(n) = mean(p(i1:i2));
end
p_mean = mean(p_mean);

for n=1:8
    t = t_all_OSSI{n};
    p = p_all_OSSI{n};

    % determine p during calibration period, and then cal coef
    i1 = find(t>cal_period_t1,1);
    i2 = find(t>cal_period_t2,1);
    cal_coef(n) = ( p_mean - mean(p(i1:i2)) ) * 1e4; % convert to Pa
end
offset_OSSI = cal_coef 

%% Plot calibration --------------
% period bucket test 1
t1a = datetime(2024,12,16,9,30,0);
t2a = datetime(2024,12,16,10,0,0);

% period bucket test 2
t1b = datetime(2024,12,23,18,0,0);
t2b = datetime(2024,12,23,20,30,0);
titles = {'bucket test 1, rbr calibrated', 'bucket test 1, rbr calibrated',...
    'bucket test 1, OSSI raw data', 'bucket test 1, OSSI raw data',...
    'bucket test 1, OSSI calibrated', 'bucket test 1, OSSI calibrated'};

figure(4), clf(4), hold on
tiledlayout(3,2,"TileSpacing","compact","Padding","compact")
for n_plot = 1:6   % Make six plots
    nexttile, hold on

    % for the upper row of plots, plot p_RBR
    if n_plot<=2
        for n=1:6
            t = t_all_RBR{n};
            p = p_all_RBR{n};
            i = 1 : 4 : round( length(t) );
            plot( t(i), p(i) )
        end

        % for the next two rows, plot raw data OSSI
    else
        calibrate = n_plot>=5 % calibrate for the last two rows;
        for n=1:8
            t = t_all_OSSI{n};
            p = p_all_OSSI{n};
            p = p + calibrate * offset_OSSI(n)/1e4;
            i = 1 : 10 : round( length(t) );

            % Plot (un)calibrated data
            if n<7,         plot( t(i), p(i)) % plot in default color
            elseif n == 7,  plot( t(i), p(i), 'c')
            else,           plot( t(i), p(i), 'm')
            end
        end
    end
    % formatting
    title( titles{n_plot} )
    grid on
    ylabel('Pressuce [ dbar]')
    if odd(n_plot)
        xlim([t1a t2a])
    else
        xlim([t1b t2b])
    end
end %end n_plot
legend

%% plot low tides/other comparison periods of OSSI sensors
% Make a plot where 12 different moments in time are plotted: bucket tests, calibration periods, and low tides (when all sensors were dry)
% NB: for comparison of different calibrations, run this section (make a plot) for every calibration perioid (first run the matlab section 
% for calibration). Mannually change subplot number. 
t0_end = [datetime(2024, 12, 16, 9,  30, 0),  datetime(2024, 12, 16, 10, 0,  0)   % bucket test 1
          datetime(2024, 12, 23, 17, 0,  0),  datetime(2024, 12, 23, 21, 0,  0)   % bucket test 2
          datetime(2024, 12, 16, 9,  43, 30), datetime(2024, 12, 16, 9,  44, 0)   % calibration period 1: 'platteau' during first bucket test
          datetime(2024, 12, 23, 19, 11, 0),  datetime(2024, 12, 23, 19, 25, 0)   % % calibration period 2: 'platteau during second bucket test
          datetime(2024, 12, 16, 20, 0,  0),  datetime(2024, 12, 17, 1,  0,  0)   % low tide 1
          datetime(2024, 12, 17, 7,  0,  0),  datetime(2024, 12, 17, 13, 0,  0)   % low tide 2
          datetime(2024, 12, 17, 20, 0,  0),  datetime(2024, 12, 18, 3,  0,  0)   % low tide 3
          datetime(2024, 12, 18, 7,  0,  0),  datetime(2024, 12, 18, 14, 0,  0)   % low tide 4
          datetime(2024, 12, 20, 10, 0,  0),  datetime(2024, 12, 20, 16, 0,  0)   % low tide 5
          datetime(2024, 12, 20, 22, 30, 0),  datetime(2024, 12, 21, 3, 30,  0)   % low tide 6
          datetime(2024, 12, 21, 11, 0,  0),  datetime(2024, 12, 21, 16, 0,  0)   % low tide 7
%           datetime(2024, 12, 23, 13, 0,  0),  datetime(2024, 12, 23, 22, 0,  0)]; % low tide 8
          datetime(2024, 12, 23, 18, 0,  0),  datetime(2024, 12, 23, 19, 0,  0)]; % cal period 3: just before bucket test 2

% Set y-lims per subplot, such that different calibrations can be plotted with the same axis, for easier comparison
ylims = [10.2 10.8
         10.2 10.8
         10.6 10.75 
         10.55 10.75
         10.15 10.4
         10.15 10.4
         10.0 10.3
         10.0 10.25
         10. 10.3
         10. 10.3
         9.95 10.25
         10.1 10.4];

figure(5), clf(5), hold on
tiledlayout(3,4,"TileSpacing",'tight','Padding','compact')
for n_plot = 1:12
    nexttile, hold on
    for n=1:8 % for every subplot, plot 8 lines (8 OSSI sensors)
        t = t_all_OSSI{n}; 
        i1 = find(t>t0_end(n_plot,1),1);
        i2 = find(t>t0_end(n_plot,2),1);

        p = p_all_OSSI{n};
        p = p(i1:10:i2) + offset_OSSI(n)/1e4;
        t = t(i1:10:i2);

        % Plot calibrated data
        if n<7,           plot( t, p)        % 6 default matlab colors
        elseif n==7,      plot( t, p, 'c')   % cyan and magenta explicitly defined
        else,             plot( t, p, 'm')
        end
    end

    grid on  % formatting of every subplot
    ylim( ylims(n_plot,:)  );
end
legend % formatting of last subplot: add legend

%% Compare calibration data
figure(6), clf(6), hold on
plot(t_all_RBR{1}, p_all_RBR{1}*100 )
plot(t_all_OSSI{8}, (p_all_OSSI{8} + offset_OSSI(8)/1e4 ) * 100 )
legend('ref.P1 RBR', 'ref.P1 OSSI')
ylabel('Atmospheric pressure [hPa]')

%% timing
% period bucket test 1
t1 = datetime(2024,12,16,9,30,0);
t2 = datetime(2024,12,16,10,0,0);

% period bucket test 1
% t1 = datetime(2024,12,16,9,49,00);  % peak/variations visible in all sensors, during bucket test
% t2 = datetime(2024,12,16,9,50,00);
% 23 dec, 13:52  % peaks accross all OSSIs 
% 23 dec, 17:52  % idem

% period bucket test 2
t1 = datetime(2024,12,23,18,0,0);
t2 = datetime(2024,12,23,20,30,0);

%%
clc
for n=1:8
    time_until_peak(n) = t_peak_OSSI(n) - t_all_OSSI{n}(1);
end
t_cal_factor = (time_until_peak + dt_OSSI) ./ time_until_peak;



%% function
function [t_select, p_select] = select_t_p(t_all,p_all, n, t_start_end)
    t_select = t_all{n}; % + dt_OSSI(n)/24/3600;
    i1 = find(t_select>t_start_end(1),1);
    i2 = find(t_select>t_start_end(2),1);

    p_select = p_all{n};
    p_select = p_select(i1:i2);
    t_select = t_select(i1:i2);
end