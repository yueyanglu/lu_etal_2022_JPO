
% calc the correlation btw two matrix along their certain dimension (e.g.,
% time)
% 

clear

homedir = getenv('HOME');
addpath(genpath([homedir '/HYCOM']));
addpath(genpath([homedir '/mytoolbox']));

%-----------------------------------------  read grid
hycom_domain = 'GSH';
read_HYCOM_grid
scp2 = scpx .* scpy;
scu2 = scux .* scuy;
scv2 = scvx .* scvy;

% mesh in [m]. Southwestern p-point is the origin [0m, 0m]
xmesh = cumsum(scux,2) - scux(:,1);
ymesh = cumsum(scvy,1) - scvy(1,:);

%% 
% layers
layers = 1:30 ; %1:15%:30%1:30;
nk = length(layers);
% times
[day_s, day_e, dt] = deal(26, 296, 5); % 26-226
t_al = day_s:dt:day_e;
nt_al = length(t_al);
d1yr = 365; % days in one year

%------------------------------- dirs
    
[UFLX_INDEX, VFLX_INDEX, DPM_INDEX, DPI_INDEX] = deal(1);
%
mld_path = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/MLD';
% mld_path = '/glade/work/yueyanglu/GSH_OUTPUT/MLD';

% 
dir_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101',...
    '/projects2/rsmas/ikamenkovich/Atlantic_HR/UVDP',...
    '/projects2/rsmas/ikamenkovich/Atlantic_HR/UVDP'};
% dir_al = {'/glade/scratch/yueyanglu/hra_expt/UVDP_CS',...
%     '/glade/scratch/yueyanglu/hra_expt/UVDP',...
%     '/glade/scratch/yueyanglu/hra_expt/UVDP'};
% /glade/scratch/yueyanglu/hra_expt/UVDP  /projects2/rsmas/ikamenkovich/Atlantic_HR/UVDP
% 

% E_al = {'/glade/scratch/yueyanglu/hra_expt/UVDP_CS/out1',...
%     '/glade/scratch/yueyanglu/hra_expt/UVDP/out3',...
%     '/glade/scratch/yueyanglu/hra_expt/UVDP/out2'};
E_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/exp0m',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/expKiso',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/expKisoA',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/expKL',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/exp2',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/exp0'};
% /glade/scratch/yueyanglu/hra_expt/UVDP/out2  /projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/run1

%
carryTracer_al = {1, 1, 1, 1, 1, 1};
ids = [1 3 4 5 6];
% ids = [1:5];
E_al = E_al(ids);  dir_al = dir_al(ids); carryTracer_al = carryTracer_al(ids);

%------------------------------- preallocate vars to be read in cell

%------- vars for vertically-averaged c 
ncel = length(E_al);
c_al = cell(1,ncel);  c_al(:) = {zeros(JDM,IDM,nk,nt_al)};

% c_al_snap = cell(1,ncel); c_al_snap(:) = {zeros(JDM,IDM,nk)};
% vertically averaged tracer
c_zsm1_al = cell(1,ncel); c_zsm1_al(:) = {zeros(JDM,IDM,nt_al)};
c_zsm2_al = cell(1,ncel); c_zsm2_al(:) = {zeros(JDM,IDM,nt_al)};
c_zsm3_al = cell(1,ncel); c_zsm3_al(:) = {zeros(JDM,IDM,nt_al)};
% c_zsm4_al = cell(1,ncel); c_zsm4_al(:) = {zeros(JDM,IDM,nt_al)};

% the # of layer whose bottom is below MLD

% func
sm_weight = @(fsum, fweight, inds, dims) squeeze( nansum(fsum.*inds, dims) ./ nansum(fweight.*inds, dims) );

%------------------------------- 
parpool('local',15);

for icel = 1:ncel %
    
    %------- what tracer experiment
    dir_path = dir_al{icel};
    E = E_al{icel};
    NTRACR = carryTracer_al{icel};
    
    fprintf('\n\nUVDP will be read from: %s ...\n', dir_path);
    fprintf('MLD will be read from: %s ...\n', mld_path);
    fprintf('Tracer-%d will be read from: %s ...\n', NTRACR, E);

    for it = 1:nt_al % [1 2 3 4]
        
        %------------------------------- current time
        nday = t_al(it);
        nyr = floor((nday - 1)/d1yr) + 1;   % start from year '0001
        yrStr = num2str(nyr, '%4.4i');
        dyStr = num2str(floor(nday - (nyr-1)*d1yr), '%3.3i');
        hrStr = num2str(mod(nday - d1yr, 1)*24, '%2.2i');
        
        fprintf('\n\nnday: %s ...\n', num2str(nday));
        
        %----------- flds of all layers at a instant
        c_alz = zeros(JDM,IDM,nk);
        cmass_alz = zeros(JDM,IDM,nk);
        
        % no NaN, only 0, same with other forcing flds . file_name
        [mld,~] = read_mld_offline_func(mld_path,dyStr,hrStr);
        fprintf(1,'Read mld...\n');

        %
        [dp_alz, uflx_alz, vflx_alz] = deal(zeros(JDM,IDM,nk)); 
        %------------------------------- all layers
        for ik = 1:nk
            % layer
            klay = layers(ik);
            fprintf('\nklay: %s ...\n', num2str(klay));

            %---- uvdp
%             read_uvdp_GSH;
            [uflx,vflx,dpm,dpi,~] = read_uvdp_GSH_func(dir_path,dyStr,hrStr,klay,UFLX_INDEX,VFLX_INDEX,DPM_INDEX,DPI_INDEX);
            fprintf('  Read UVDP.\n');
            fland = dpi <= 1.e-10;
            dpi(fland) = NaN;
            dp_alz(:,:,ik) = dpi;
            uflx_alz(:,:,ik) = uflx;
            vflx_alz(:,:,ik) = vflx;
            %---- c
%             diag_tracer;
            [tracer,~] = diag_tracer_func(E,yrStr,dyStr,hrStr,klay,NTRACR);
            fprintf('  Read tracer-%d. %s\n',NTRACR);
            tracer(fland) = NaN;
            tracer(tracer < 1e-9) = NaN;
            
            %----- c/cmass of all layers at current instant
            c_alz(:,:,ik) = tracer;
            cmass_alz(:,:,ik) = tracer .* dpi .* scp2;
%             %----- total c-mass of each layer
%             cmass_al{icel}(ik,it) = nansum(tracer .* dpi .* scpx.*scpy, 'all');
            %----- 
%             c_al{icel}(:,:,ik,it) = tracer;

        end % nk
        
        % Works only if nk = KDM
        if nk == KDM
            %----- depth of each layer bottom [j-i-k]
            depth_alz = cumsum(dp_alz, 3);
            %----- if the layer' bottom is below MLD
            ifAboveMLD_3d = depth_alz <= repmat(mld, [1 1 KDM]); % [j,i,k]
            ifAboveMLD_3d(:,:,1) = 1; % make sure at least one layer is counted!
            % above middle layer
            k2 = 15;
            k1d = [ones(1,k2), zeros(1,KDM-k2)];
            ifAboveMidle_3d = repmat(reshape(k1d,[1 1 KDM]), [JDM IDM 1]);

            %----- if middle layer: below MLD ~ above k2 (included)
            ifSec2_3d = ~ifAboveMLD_3d .* ifAboveMidle_3d;
            %----- if bottom layer: below MLD && below k2 (< k2) 
            ifSec3_3d = ~ifAboveMidle_3d .* ~ifAboveMidle_3d;
            
            %----- vertically averaged c in different section [j,i,t]
            c_zsm1_al{icel}(:,:,it) = sm_weight(cmass_alz, dp_alz .* scp2, ifAboveMLD_3d, [3]);
            c_zsm2_al{icel}(:,:,it) = sm_weight(cmass_alz, dp_alz .* scp2, ifSec2_3d, [3]);
            c_zsm3_al{icel}(:,:,it) = sm_weight(cmass_alz, dp_alz .* scp2, ifSec3_3d, [3]);
%             c_zsm4_al{icel}(:,:,it) = sm_weight(cmass_alz, dp_alz .* scp2, ifAboveMidle_3d, [3]);

        end

    end %icel
end
delete(gcp('nocreate'))

% 'c_al',　'\<c_zsm\w*'
save(['c_evo_zsm2.mat'],'-regexp','\<c_zsm\w*','E_al','t_al','nt_al','nk','layers','-v7.3');

%% stats
%---- 
dim = 3;
% corr matrix at each layer, M&F, ADV&F
[corr13_alz, corr23_alz] = deal(NaN*zeros(JDM,IDM,KDM));
% 
% [corr13_zsm1, corr23_zsm1, corr13_zsm3, corr23_zsm3] = deal(NaN*zeros(JDM,IDM));

% Frobenius norm of the tracer differece matrix at each layer and each time
[normF13, normF23] = deal(NaN*zeros(KDM,nt_al));

for ik = 1:KDM
    disp(['ik = ' num2str(ik)]);
    % J-I-T
    f1 = squeeze(c_al{1}(:,:,ik,:));
    f2 = squeeze(c_al{2}(:,:,ik,:));
    f3 = squeeze(c_al{3}(:,:,ik,:));
    
    %--- corr matrix
    corr13_alz(:,:,ik) = corr_array(f1,f3,dim);
    corr23_alz(:,:,ik) = corr_array(f2,f3,dim);
    
    %--- Frobenius norm
    for it = 1:nt_al
        dist13 = f1(:,:,it) - f3(:,:,it);
        dist23 = f2(:,:,it) - f3(:,:,it);
        % nan
        dist13(isnan(dist13)) = [];
        dist23(isnan(dist23)) = [];
        % norm
        normF13(ik,it) = norm(dist13,'fro');
        normF23(ik,it) = norm(dist23,'fro');
    end
end

% for vertically averaged tracer
corr13_zsm1 = corr_array(c_zsm1_al{1},c_zsm1_al{3},dim);
corr23_zsm1 = corr_array(c_zsm1_al{2},c_zsm1_al{3},dim);
% 
corr13_zsm3 = corr_array(c_zsm3_al{1},c_zsm3_al{3},dim);
corr23_zsm3 = corr_array(c_zsm3_al{2},c_zsm3_al{3},dim);

save('stats_cEvo_2','corr13_*','corr23_*','normF13','normF23','E_al','t_al','nt_al');

%% stats 2
 
dim = 3;
icel_test = [1 2 3 4]; nc_test = length(icel_test);
icel_bench = 5;

%-------- Frobenius norm of the tracer differece matrix vs time
normF_zsm1 = cell(1,nc_test);
normF_zsm3 = cell(1,nc_test);
% 
normF_zsm1(:) = {zeros(1,nt_al)};
normF_zsm3(:) = {zeros(1,nt_al)};
% 
for icel = 1:nc_test
    icel_do = icel_test(icel);
    %--- Frobenius norm at each time
    for it = 1:nt_al
        diff_zsm1 = c_zsm1_al{icel_do}(:,:,it) - c_zsm1_al{icel_bench}(:,:,it);
        diff_zsm3 = c_zsm3_al{icel_do}(:,:,it) - c_zsm3_al{icel_bench}(:,:,it);
        % nan
        diff_zsm1(isnan(diff_zsm1)) = [];
        diff_zsm3(isnan(diff_zsm3)) = [];
        % norm
        normF_zsm1{icel}(1,it) = norm(diff_zsm1,'fro');
        normF_zsm3{icel}(1,it) = norm(diff_zsm3,'fro');
    end
end

%-------- 2D corr matrix for vertically averaged tracer

corr_zsm1 = cell(1,nc_test);
corr_zsm3 = cell(1,nc_test);
% 
for icel = 1:nc_test
    icel_do = icel_test(icel);
    corr_zsm1{icel} = corr_array(c_zsm1_al{icel_do},c_zsm1_al{icel_bench},dim);
    corr_zsm3{icel} = corr_array(c_zsm3_al{icel_do},c_zsm3_al{icel_bench},dim);
end

save('stats_cEvo_2','corr_zsm*','normF_zsm*','E_al','t_al','nt_al');

%% plot

% load('stats_cEvo.mat')

% ----------- corr map
ik = 15;
icel = 4;
titlstr = {'$\textrm{Corr}( c_{MEAN}, c_{FULL} )$', '$\textrm{Corr}( c_{ADV}, c_{FULL} )$'};
f_plt = {corr_zsm1{icel}, corr_zsm3{icel}};

font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
[ha, ~] = tight_subplot(1,2,[.05 .03],[.08 .05],[.08 .08]);
% 
for icel = 1:2
    axes(ha(icel));
    plot_field_model(f_plt{icel},plon1d,plat1d,'balance');
    cmocean('balance',10)
    caxis([-1 1])
    if icel == 1
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,...
            'ytick',30:5:45,'linewidth',1.2);
    else
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,'xticklabels',[],...
            'ytick',30:5:45,'yticklabels',[],'linewidth',1.2);
    end
    ax = gca;
    ax.Title.Interpreter = 'Latex';
    ax.Title.String = titlstr{icel};
    ax.Title.FontSize = 16;
end
% 
cb = colorbar;
set(cb,'orientation','horizontal','Position',[0.35 0.26 0.30 0.02]);

%% ----------- norm of the trracer difference, z vs. time

clim = [0 140];
it_plt = 1:6:nt_al;
cmap = 'amp';
f_plt = {normF13, normF23};
% f_plt = {normF_zsm1{1}, normF_zsm1{2}};
% titlstr = {'$|\!|\textrm{MEAN} - \textrm{FULL}|\!|_F$', '$|\!|\textrm{EXP-ADV} - \textrm{FULL}|\!|_F$'};
titlstr = {'$|\!| c_{MEAN} - c_{FULL} |\!|_F$', '$|\!| c_{ADV} - c_{FULL} |\!|_F$'};

font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
[ha, ~] = tight_subplot(2,2,[.10 .04],[.12 .08],[.10 .10]);
% 
for icel = 1:2
    axes(ha(icel));
%     [~,h] = contourf(1:nt_al, 1:KDM, f_plt{icel}, 16); set(h,'EdgeColor','none');
    h = pcolor(1:nt_al, 1:KDM, f_plt{icel}); set(h,'EdgeColor','none');
    caxis(clim);
    cmocean(cmap)
    ax = gca;
    ax.XLim = [1 nt_al];
    ax.XTick = it_plt;  ax.XTickLabel = t_al(it_plt) - t_al(1) ;
    ax.YLim = [1 KDM];
    ax.YDir = 'reverse';
    ax.TickLength = [.01, .01];  ax.LineWidth = 1.5; ax.TickDir = 'out';
    ax.XAxis.FontSize = 12; ax.YAxis.FontSize = 12;
    % label
    ax.XAxis.Label.String = 'Day';
    ax.YAxis.Label.String = 'Layer';
    % title
    ax.Title.FontSize = 16;
    ax.Title.Interpreter = 'latex';
    ax.Title.String = titlstr{icel};
    %
    if icel == 2 
        ax.YTickLabel = '';
        ax.YAxis.Label.String = '';
    end
end
% 
cb = colorbar;
set(cb,'Position',[0.92 0.33 0.02 0.40]);

%%
figure
for icel = 1:4
    plot(normF_zsm3{icel})
    hold on
end
legend('MEAN', 'KisoA', 'KL', 'EXP-ADV')

%% plot Frob norm + 4 corr matrix in mixed layer

x_top = 0.105; y_top = 0.72; h_top = 0.26; w_top = 0.67; 
x_bot = 0.08; y_bot = 0.04; h_bot = 0.27; w_bot = 0.36; dx_bot = 0.001; dy_bot = 0.041;

pos_top = [x_top y_top w_top h_top];

% -------------- plt top panel (Frob norm)
f_plt_top = normF_zsm1;
ncel = 4;
ylabel_norm = '$|\!| c - c_{FULL} |\!|_F$';
titlstr = {'$\textrm{MEAN}$', '$\textrm{EXP-}\textbf{K}_{red}$', ...
    '$\textrm{EXP-}\kappa \mbox{\boldmath $\chi$} $', '$\textrm{EXP-ADV}$'}; % EXP-ADV
linewidths = 2;
linetyles = {'k--','b-','r-','g-'};
firstDay = 26;
relative_times = t_al - firstDay + 1; % 
its = [1:5:nt_al nt_al];

% -- plt 
font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
% subplot('Position',pos_top);
for icel = 1:ncel
    subplot('Position',pos_top);
    plot(1:nt_al,f_plt_top{icel}, linetyles{icel},'LineWidth',linewidths)
    hold on
end
ax = gca;
ax.TickLength = [.01, .01];
ax.LineWidth = 1.0;
ax.XLim = [0 nt_al+1];
ax.XTick = its;
ax.XTickLabel = relative_times(its);
ax.XAxis.Label.String = 'Day';
ax.YAxis.Label.String = ylabel_norm; ax.YAxis.Label.Interpreter = 'Latex'; 
ax.YAxis.Label.FontSize = 10;
lgd = legend(titlstr,'Interpreter','latex','FontSize',6);
lgd.Location = 'northeast'; lgd.Box = 'on';
hold off

% -------------- plt bottom panel (corr mat 2X2)
f_plt_bot = corr_zsm1;
titlstr_bot = {'$\textrm{Corr}( \textrm{MEAN}, \textrm{FULL} )$', ...
    '$\textrm{Corr}(\textrm{EXP-}\textbf{K}_{red}, \textrm{FULL})$',...
   '$\textrm{Corr} (\textrm{EXP-}\kappa \mbox{\boldmath $\chi$}, \textrm{FULL})$',...
   '$\textrm{Corr}( \textrm{EXP-ADV}, \textrm{FULL} )$'};

pos_bot = { [x_bot y_bot+h_bot+dy_bot w_bot h_bot], [x_bot+w_bot+dx_bot y_bot+h_bot+dy_bot w_bot h_bot], ...
    [x_bot y_bot w_bot h_bot], [x_bot+w_bot+dx_bot y_bot w_bot h_bot] };

% -- plt 
for icel = 1:ncel
    subplot('Position',pos_bot{icel});
%     subplot(2,2,icel);
    plot_field_model(f_plt_bot{icel},plon1d,plat1d,'balance');
    cmocean('balance',10)
    caxis([-1 1])
    if icel == 3
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,...
            'ytick',30:5:45,'linewidth',1.2);
    else
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,'xticklabels',[],...
            'ytick',30:5:45,'yticklabels',[],'linewidth',1.2);
    end
    ax = gca;
    ax.Title.Interpreter = 'Latex';
    ax.Title.String = titlstr_bot{icel};
    ax.Title.FontSize = 10;
end

cb = colorbar;
set(cb,'Position',[x_bot+2*w_bot-0.01 y_bot 0.02 h_bot]);

% -------- plot (a), (b)...
dim1 = [0.11 0.91 0.06 0.06];
dim2 = [0.11 0.55 0.06 0.06];
dim3 = [0.11+w_bot+dx_bot 0.55 0.06 0.06];
dim4 = [0.11 0.24 0.06 0.06];
dim5 = [0.11+w_bot+dx_bot 0.24 0.06 0.06];

annotation('textbox',dim1,'String','(a)','FitBoxToText','on','EdgeColor','none','BackgroundColor','none','fontsize',14);
annotation('textbox',dim2,'String','(b)','FitBoxToText','on','EdgeColor','none','BackgroundColor','none','fontsize',14);
annotation('textbox',dim3,'String','(c)','FitBoxToText','on','EdgeColor','none','BackgroundColor','none','fontsize',14);
annotation('textbox',dim4,'String','(d)','FitBoxToText','on','EdgeColor','none','BackgroundColor','none','fontsize',14);
annotation('textbox',dim5,'String','(e)','FitBoxToText','on','EdgeColor','none','BackgroundColor','none','fontsize',14);

% printpdf(gcf,'fig10','-r600')

%% --- two metrics together

% --- plt top panel
clim = [0 140];
it_plt = 1:12:nt_al;
cmap = 'amp';
f_plt_top = {normF13, normF23};
titlstr_top = {'$|\!| c_{MEAN} - c_{FULL} |\!|_F$', '$|\!| c_{ADV} - c_{FULL} |\!|_F$'};

% --- plt bottom panel
ik = 15;
titlstr_bot = {'$\textrm{Corr}( c_{MEAN}, c_{FULL} )$', '$\textrm{Corr}( c_{ADV}, c_{FULL} )$'};
f_plt_bot = {corr13_alz(:,:,ik), corr23_alz(:,:,ik)};

% ---- positionns
x_top = 0.13; y_top = 0.5; h_top = 0.44; w_top = 0.3; dx_top = 0.02;
x_bot = 0.08; y_bot = 0.04; h_bot = 0.30; w_bot = 0.35; dx_bot = 0.02;

pos_top = { [x_top y_top w_top h_top], [x_top+w_top+dx_top y_top w_top h_top] };
pos_bot = { [x_bot y_bot w_bot h_bot], [x_bot+w_bot+dx_bot y_bot w_bot h_bot] };


font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
% --- top
for icel = 1:2
    subplot('Position',pos_top{icel});
%     [~,h] = contourf(1:nt_al, 1:KDM, f_plt{icel}, 16); set(h,'EdgeColor','none');
    h = pcolor(1:nt_al, 1:KDM, f_plt{icel}); set(h,'EdgeColor','none');
    caxis(clim);
    cmocean(cmap)
    ax = gca;
    ax.XLim = [1 nt_al];
    ax.XTick = it_plt;  ax.XTickLabel = t_al(it_plt) - t_al(1) ;
    ax.YLim = [1 KDM];
    ax.YDir = 'reverse';
    ax.TickLength = [.01, .01];  ax.LineWidth = 1.5; ax.TickDir = 'out';
    ax.XAxis.FontSize = 12; ax.YAxis.FontSize = 12;
    % label
    ax.XAxis.Label.String = 'Day';
    ax.YAxis.Label.String = 'Layer';
    % title
    ax.Title.FontSize = 16;
    ax.Title.Interpreter = 'latex';
    ax.Title.String = titlstr_top{icel};
    %
    if icel == 2 
        ax.YTickLabel = '';
        ax.YAxis.Label.String = '';
    end
end
% 
cb = colorbar;
set(cb,'Position',[.76 y_top 0.02 h_top]);

% --- bot
for icel = 1:2
    subplot('Position',pos_bot{icel});
    plot_field_model(f_plt_bot{icel},plon1d,plat1d,'balance');
    cmocean('balance',10)
    caxis([-1 1])
    if icel == 1
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,...
            'ytick',30:5:45,'linewidth',1.2);
    else
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,'xticklabels',[],...
            'ytick',30:5:45,'yticklabels',[],'linewidth',1.2);
    end
    ax = gca;
    ax.Title.Interpreter = 'Latex';
    ax.Title.String = titlstr_bot{icel};
    ax.Title.FontSize = 16;
end
% 
cb = colorbar;
set(cb,'Position',[.8 y_bot 0.02 h_bot]);

% -------- plot (a), (b)...
dim1 = [0.13 0.50 0.06 0.06];
dim2 = [0.13+w_top+dx_top 0.50 0.06 0.06];
dim3 = [0.1 0.27 0.06 0.06];
dim4 = [0.1+w_bot+dx_bot 0.27 0.06 0.06];

annotation('textbox',dim1,'String','(a)','FitBoxToText','on','EdgeColor','none','BackgroundColor','none','fontsize',14);
annotation('textbox',dim2,'String','(b)','FitBoxToText','on','EdgeColor','none','BackgroundColor','none','fontsize',14);
annotation('textbox',dim3,'String','(c)','FitBoxToText','on','EdgeColor','none','BackgroundColor','none','fontsize',14);
annotation('textbox',dim4,'String','(d)','FitBoxToText','on','EdgeColor','none','BackgroundColor','none','fontsize',14);


% set(gcf,'PaperPositionMode','auto'); print(gcf,'fig_metrcs_3exp','-dpng','-r600');
