clear, clc

map = 'O:\HybridDune experiment\data RBR, OSSI\copy RBR Udrive series1\QC\';
naam = {'S2P3 RBR1 p_rel.nc'
        'S4P3 RBR2 p_rel.nc'
        'S1P2 RBR3 p_rel.nc'
        %'refP1 RBR4 p_rel.nc'
        'S1P3 RBR5 p_rel.nc'
        'S3P3 RBR6 p_rel.nc'};
t_factor = [1e3 1e3 1e3 1e6 1e6];
rho = 1027;
g = 9.81;


for n_file = 3
    % Load data
    p0 = ncread([map,naam{n_file}],'p_rel');
    p0_full = reshape(p0,[],1);
    t = ncread([map,naam{n_file}],'t'); % t per 20 minutes
    t1 = datetime(2024,12,12,9,0,0) + double(t)/t_factor(n_file)/3600/24; % Make datetime from (mili)seconds
    t = ncread([map,naam{n_file}],'N');  % t within block. In seconds
    t_2d = t1' + t/3600/24;  % 1 column per block
    t_full = reshape(t_2d,[],1); % 1 column for entire timeseries
    
    %% filter Ruben
    clc
    p0_ruben = p0;
    p0_ruben(p0_ruben < 0.3*rho*g) = nan;
    sum(std(p0_ruben, 'omitnan') < 70);  % count number of blocks where applicable
    p0_ruben(:,std(p0_ruben, 'omitnan') < 70) = nan;

    p0_ruben2 = p0_ruben;
    sum(abs(p0_ruben2 - mean(p0_ruben2,'omitnan')) > 3 * std(p0_ruben2,'omitnan')); % count number of obs where applicable
    p0_ruben2( abs(p0_ruben2 - mean(p0_ruben2,'omitnan')) > 3 * std(p0_ruben2,'omitnan') ) = nan;
    p0_ruben2_full = reshape(p0_ruben2,[],1);

    % filter Daan
    p0_daan = p0;
    p0_daan(:, sum(p0_daan<0.05*rho*g) > 0.1 * size(p0_daan,1)) = nan;
    p0_daan2 = p0_daan;
    p0_daan2(p0_daan2<0)=0;
    p0_daan2_full = reshape(p0_daan2,[],1);

    p0_daan3 = p0;
    p0_daan3(:, sum(p0_daan3<0.05*rho*g) > 0.25 * size(p0_daan3,1)) = nan;
    p0_daan3_full = reshape(p0_daan3,[],1);
    %%
    clc
    y0 = p0(:,599)/1e4;
    y = p0_ruben(:,599)/1e4;
    y2 = p0_daan(:,599)/1e4;

    mean(y,'omitnan')
    std(y,'omitnan')
    mean(y,'omitnan') - 3 *     std(y,'omitnan')
    mean(y,'omitnan') + 3 *     std(y,'omitnan')

    figure(1), clf(1),
    subplot(1,2,1), hold on
    plot(y0)
    plot(y)
    plot([0 1e4],    [1 1]*mean(y0,'omitnan') - 3 *     std(y0,'omitnan'))
    plot([0 1e4],    [1 1]*mean(y0,'omitnan') + 3 *     std(y0,'omitnan'))

    subplot(1,2,2), hold on
    plot(y0)
    plot(y2)
    [mean(y0),std(y0)
     mean(y, 'omitnan'),std(y,'omitnan')
     mean(y2),std(y2)]
    %% make plot
    figure(2), clf(2), hold on
%     tiledlayout(2,2)

    t_start1 = datetime(2024,12,17,10,0,0);
    t_end1   = datetime(2024,12,23,16,0,0);
    t_start2 = datetime(2024,12,20,6,0,0);
    t_end2   = datetime(2024,12,20,7,0,0);

    i_start1 = find(t_full>t_start1,1);
    i_end1   = find(t_full>t_end1,1);
    i_start2 = find(t_full>t_start2,1);
    i_end2   = find(t_full>t_end2,1);

%     nexttile, hold on
%     plot(t_full(i_start1:i_end1),p0_full(i_start1:i_end1)/1e4, '.')
%     plot(t_full(i_start1:i_end1),p0_ruben2_full(i_start1:i_end1)/1e4, '.')
%     plot(xlim,[0 0],'k')
%     ylabel('h above sensor [m]')
%     title('filter Ruben')
% 
%     nexttile, hold on
%     plot(t_full(i_start2:i_end2),p0_full(i_start2:i_end2)/1e4, '.')
%     plot(t_full(i_start2:i_end2),p0_ruben2_full(i_start2:i_end2)/1e4, '.')
%     plot(xlim,[0 0],'k')
%     ylabel('h above sensor [m]')
%     title('filter Ruben')

%     nexttile, hold on
    plot(t_full(i_start1:i_end1),p0_full(i_start1:i_end1)/1e4, '.')
    plot(t_full(i_start1:i_end1),p0_daan3_full(i_start1:i_end1)/1e4, '.')
    plot(t_full(i_start1:i_end1),p0_daan2_full(i_start1:i_end1)/1e4, '.')
    plot(xlim,[0 0],'k')
%     title('filter Daan')

%     nexttile, hold on
%     plot(t_full(i_start2:i_end2),p0_full(i_start2:i_end2)/1e4, '-')
%     plot(t_full(i_start2:i_end2),p0_daan2_full(i_start2:i_end2)/1e4, '.')
%     plot(xlim,[0 0],'k')
%     title('filter Daan')

end
