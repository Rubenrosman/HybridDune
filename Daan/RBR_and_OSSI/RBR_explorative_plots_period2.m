%% import RBR NetCDF data
clear, clc

map = 'O:\HybridDune experiment\data RBR, OSSI\copy RBR Udrive series2\raw NetCDF\';
naam = {'refP1 RBR3.nc'
        'S3P3 RBR6.nc'
        'S1P2 RBR2.nc'
        'S3P2 RBR4.nc'
        'S4P2 RBR1.nc'
        'S3P1 RBR5.nc'};
F_all = [8 16 8 8 8 16];
t_all = [];
p_all = [];

for n=1:6
    p = ncread([map,naam{n}],'p')/1e5;
    t = ncread([map,naam{n}],'t'); % time 
    t_atts = ncinfo([map,naam{n}],'t').Attributes; % unit of time
    t = read_time_from_xarray_netcdf(t,t_atts);
    
    t_all{n} = datetime(2024,12,12,9,0,0) + t_days;
    p_all{n} = p;
    legend_labels{n}=naam{n}(1:end-3);
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
figure(3), clf(3),hold on
plot(t_all{1}(i1:i2),(p_matrix-0*delta_p_cal)*1000)
title('Calibration period')
ylabel('Pressure [mbar = cm water equiv.]')
legend

% Plot calibrated data
titles = {'Calibrated, before employment','Calibrated, after employment','Calibrated, bucket test 1','Calibrated, bucket test 2'};
xlims = [datetime(2024,12,12,0, 0, 0), datetime(2024,12,17,0, 0,0);
         datetime(2024,12,23,18,30,0), datetime(2024,12,27,0, 0,0);
         datetime(2024,12,16,8, 30,0), datetime(2024,12,16,11,0,0);
         datetime(2024,12,23,18,30,0), datetime(2024,12,23,21,0,0)];
ylims=[1025 1050;
       1030 1055;
       1037 1082;
       1030 1075]

fig=figure(4); clf(4)
% fig.Position = [100 100 800 700];
tiledlayout(2,2,'TileSpacing','compact','Padding','tight')
p_all_filtered = p_all;
p_all_filtered{3}( (t_all{3}>datetime(2024,12,16,13,0,0)) & (t_all{3}<datetime(2024,12,16,20,0,0)) ) = nan;
for i=1:4
    ax=nexttile; hold on
    title(titles(i))
    for n=1:6
        plot(t_all{n}(1:5:end), (p_all_filtered{n}(1:5:end)-delta_p_cal(n))*1000 )
    end
    xlim(xlims(i,:))
    ylabel('Pressure [mbar = cm water equiv.]')
    ylim(ylims(i,:))
end
L=legend(legend_labels,'Location','NEC');

%% Plot reference data
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

%% Plot reference data plus sensor 3
fig=figure(5); clf(5),hold on
tiledlayout(2,2, 'Padding','tight')
xlims=[datetime(2024,12,16,20,0,0), datetime(2024,12,16,21,0,0);
       datetime(2024,12,16,20,0,0), datetime(2024,12,16,20,30,0);
       datetime(2024,12,16,20,0,0), datetime(2024,12,16,20,5,0)];
window_seconds = 100;
for n=1:3
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
fig=figure(5); clf(5);  hold on
% fig.Position = [100 100 800 700];
p_all_filtered = p_all;
p_all_filtered{3}( (t_all{3}>datetime(2024,12,16,13,0,0)) & (t_all{3}<datetime(2024,12,16,20,0,0)) ) = nan;

for n=1:6
    y = movmean(p_all{n} - delta_p_cal(n), 60*8 ) * 1000;
    plot(t_all{n}(1:60:end), y(1:60:end) )
end
ylabel('Pressure [mbar = cm water equiv.]')

L=legend(legend_labels,'Location','NEC');