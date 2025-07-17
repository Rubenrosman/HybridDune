%% import RBR NetCDF data
clear, clc

map = 'O:\HybridDune experiment\data RBR, OSSI\copy RBR Udrive series1\raw NetCDF\';
naam = {'S2P3 RBR1 raw data - period 1.nc'
        'S4P3 RBR2 raw data - period 1.nc'
        'S1P2 RBR3 raw data - period 1.nc'
        'refP1 RBR4 raw data - period 1.nc'
        'S1P3 RBR5 raw data - period 1.nc'
        'S3P3 RBR6 raw data - period 1.nc'};
F_all = [8 8 8 8 16 16];
t_all = [];
p_all = [];
for n=1:6
    p = ncread([map,naam{n}],'p')/1e5;
    t = ncread([map,naam{n}], 't');                    % load time + attributes of time
    t_atts = ncinfo([map,naam{n}],'t').Attributes;      % Load attributes of time, select unit used
    t = read_time_from_xarray_netcdf(t,t_atts);         % convert to datetime format

    t_all{n} = t;
    p_all{n} = p;
    legend_labels{n}=naam{n}(1:end-23);
end
% a = ncinfo([map,naam{1}])
%%
for n=1:6
    t1(n,1) = t_all{n}(1);
    t_end(n,1) = t_all{n}(end);
end
%% Plot all six lines
fig=figure(1); clf(1),hold on
fig.Position = [100 100 800 500];
for n=1:6
    plot(t_all{n}, p_all{n}*1000)
end
title('All six RBR sensors, raw data')
ylabel('Pressure [mbar = cm water equiv.]')
legend(legend_labels)
ax=gca;
ax.Position(1) = 0.07;
ax.Position(3) = 0.90;

%% plot reference sensor, bucket tests
fig=figure(2); clf(2)
fig.Position = [100 100 800 500];
tiledlayout(2,2,'TileSpacing','compact','Padding','tight')
nexttile(1,[1 2])
plot(t_all{4}, p_all{4}*1000)
ylabel('Pressure [mbar = cm water equiv.]')
title('Reference sensor')

nexttile,hold on
title('Bucket test 1, all sensors')
for n=1:6
    plot(t_all{n}, (p_all{n})*1000 )
end
xlim([datetime(2024,12,16,8,30,0), datetime(2024,12,16,11,0,0)])
ylabel('Pressure [mbar = cm water equiv.]')
L=legend(legend_labels,'Location','NEC');

nexttile,hold on
title('Bucket test 2, all sensors')
for n=1:6
    plot(t_all{n}, (p_all{n})*1000 )
end
xlim([datetime(2024,12,23,18,30,0), datetime(2024,12,23,21,0,0)])
ylabel('Pressure [mbar = cm water equiv.]')
% legend(legend_labels)

%% calibrate
clc
p_matrix=[]
i1 = find( t_all{1} > datetime(2024,12,12,20,0,0), 1);
i2 = find( t_all{1} > datetime(2024,12,13,6,0,0), 1);
i1 = find( t_all{1} > datetime(2024,12,23,18,0,0), 1);      % just before bucket test 2. New period: also possible for OSSIs
i2 = find( t_all{1} > datetime(2024,12,23,19,0,0), 1);
i1 = find( t_all{1} > datetime(2024,12,23,19,11,0), 1);  % 'platteau' during second bucket test. Also possible for OSSIs. Good performance for both
i2 = find( t_all{1} > datetime(2024,12,23,19,25,0), 1);
% i1 = find( t_all{1} > datetime(2024,12,15,0,0,0), 1);   % worse performance, especially for the reference sensor
% i2 = find( t_all{1} > datetime(2024,12,16,0,0,0), 1);
% 
% i1 = find( t_all{1} > datetime(2024,12,23,20,3,0), 1);  % similar performance. But short period (30min) and in (mostly drained) bucket, so
% i2 = find( t_all{1} > datetime(2024,12,23,20,33,0), 1); % not true air pressure
for n=1:6
    if n<5
        p_matrix(:,n) = [p_all{n}(i1:i2)];
    else
        p_matrix(:,n) = [p_all{n}(i1*2-1 : 2 : i2*2-1)];
    end
end
delta_p_cal = mean(p_matrix) - mean(p_matrix(:));
[1:6; -round(delta_p_cal*1e5)]


% Plot calibrated data
titles = {'Raw, before installation at beach','Raw, after removal from beach','Raw, bucket test 1','Raw, bucket test 2';
   'Calibrated, before installation at beach','Calibrated, after removal from beach','Calibrated, bucket test 1','Calibrated, bucket test 2'};

xlims = [datetime(2024,12,12,0, 0, 0), datetime(2024,12,17,0, 0,0);
         datetime(2024,12,23,18,30,0), datetime(2024,12,27,0, 0,0);
         datetime(2024,12,16,8, 30,0), datetime(2024,12,16,11,0,0);
         datetime(2024,12,23,18,30,0), datetime(2024,12,23,21,0,0)];
ylims=[1020 1050;
       1020 1050;
       1025 1080;
       1020 1075];

for n_fig = 4
fig=figure(15); clf(15)
fig.Position = [100 100 800 700];
tiledlayout(2,2,'TileSpacing','compact','Padding','tight')
p_all_filtered = p_all;
p_all_filtered{3}( (t_all{3}>datetime(2024,12,16,13,0,0)) & (t_all{3}<datetime(2024,12,16,20,0,0)) ) = nan;
for i=1:4
    ax=nexttile; hold on
    title(titles(n_fig-2,i))
    for n=1:6
        if n_fig == 3
            plot(t_all{n}(1:5:end), (p_all_filtered{n}(1:5:end))*1000 )
        else
            plot(t_all{n}(1:5:end), (p_all_filtered{n}(1:5:end)-delta_p_cal(n))*1000 )
        end
    end
    xlim(xlims(i,:))
    ylabel('Pressure [mbar = cm water equiv.]')
    ylim(ylims(i,:))
end
L=legend(legend_labels,'Location','NEC');
end

%% Plot air pressure data
fig=figure(5); clf(5),hold on
tiledlayout(1,1, 'Padding','tight')
nexttile, hold on
fig.Position = [100 100 800 500];
p = p_all{4}*1000;
plot(t_all{4}, p)
plot(t_all{4}, movmean(p,80))

title('RBR reference sensor')
ylabel('Pressure [mbar = cm water equiv.]')
legend('Raw data','10s moving avg')
xlim([datetime(2024,12,20,5,0,0), datetime(2024,12,20,12,0,0)])

%% Plot reference data plus sensor 3, to check timescale of changes and correlation between sensors
fig=figure(6); clf(6),hold on
tiledlayout(2,2, 'Padding','tight')
xlims=[datetime(2024,12,16,20,0,0), datetime(2024,12,16,21,0,0);
       datetime(2024,12,16,20,0,0), datetime(2024,12,16,20,30,0);
       datetime(2024,12,16,20,0,0), datetime(2024,12,16,20,5,0)];
window_seconds = 100;
for n=1:3 % 3 subplots, zooming in on increasingly short periods
nexttile, hold on
fig.Position = [100 100 800 500];
p = ( p_all{4} - delta_p_cal(4) ) * 1000 + 0.5;
plot(t_all{4}, movmean(p,window_seconds*10))

p = ( p_all{3} - delta_p_cal(3) ) * 1000;
plot(t_all{3}, movmean(p,window_seconds*10))

title('RBR reference sensor')
ylabel('Pressure [mbar = cm water equiv.]')
legend('RBR S1P2','RBR ref')
% xlim([datetime(2024,12,16,20,30,0), datetime(2024,12,16,21,0,0)])
xlim(xlims(n,:))
end

%% plot all sensors, calibrated, moving average
%% Plot calibrated data
fig=figure(7); clf(7);  hold on
% fig.Position = [100 100 800 700];
p_all_filtered = p_all;
p_all_filtered{3}( (t_all{3}>datetime(2024,12,16,13,0,0)) & (t_all{3}<datetime(2024,12,16,20,0,0)) ) = nan;

for n=1:6
    y = movmean(p_all{n} - delta_p_cal(n), 60*8 ) * 1000;
    plot(t_all{n}(1:60:end), y(1:60:end) )
end
ylabel('Pressure [mbar = cm water equiv.]')

L=legend(legend_labels,'Location','NEC');