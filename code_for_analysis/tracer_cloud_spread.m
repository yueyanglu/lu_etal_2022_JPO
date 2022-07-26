
addpath(genpath('/nethome/yxl1496/HYCOM'));
addpath(genpath('/nethome/yxl1496/mytoolbox'));

clear

%-----------------------------------------  read grid
hycom_domain = 'GSH';
read_HYCOM_grid
scp2 = scpx .* scpy;
scu2 = scux .* scuy;
scv2 = scvx .* scvy;

%% 
% layers
layers =  1:30;
nk = length(layers);
% times
[day_s, day_e, dt_save] = deal(22, 316, 5); % 26-271
t_al = day_s:dt_save:day_e;
nt_al = length(t_al);
d1yr = 365; 

%------------------------------- dirs
    
[UFLX_INDEX, VFLX_INDEX, DPM_INDEX, DPI_INDEX] = deal(1);
%
dir_al = {'/scratch/projects/ome/hra_expt/UVDP_CS',... 
    '/scratch/projects/ome/hra_expt/UVDP_CS',...
    '/scratch/projects/ome/hra_expt/UVDP',...
    '/scratch/projects/ome/hra_expt/UVDP'};
% 
carryTracer_al = {1, 1, 1, 1};
E_al = {'/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/exp0m',...
    '/scratch/projects/ome/hra_expt/UVDP_CS/out5',...
    '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/exp0',...
    '/scratch/projects/ome/hra_expt/UVDP/out1'};
ids = [1 2 3 ];
dir_al = dir_al(ids);
E_al = E_al(ids);

%------------------------------- preallocate vars to be read in cell 
ncel = length(E_al);
%-------
[intc_al, varxx_al, varxy_al, varyy_al, lmd1_al, lmd2_al, phi_al] = deal(cell(1,ncel));
% c_al = cell(1,ncel);
% c_al(:) = {zeros(JDM,IDM,nk,nt_al)};
% assign
intc_al(:) = {zeros(nk,nt_al)};
varxx_al(:) = {zeros(nk,nt_al)};
varxy_al(:) = {zeros(nk,nt_al)};
varyy_al(:) = {zeros(nk,nt_al)};
% 
lmd1_al(:) = {zeros(nk,nt_al)};
lmd2_al(:) = {zeros(nk,nt_al)};
phi_al(:) = {zeros(nk,nt_al)};

%------------------------------- each time
% mesh in [m]. Southwestern p-point is the origin [0m, 0m]
xmesh = cumsum(scux,2) - scux(:,1);
ymesh = cumsum(scvy,1) - scvy(1,:);

for it = 1:nt_al
    
    %------------------------------- current time
    nday = t_al(it);
    nyr = floor((nday - 1)/d1yr) + 1;   % start from year '0001
    year_num = num2str(nyr, '%4.4i'); 
    day_num = num2str(floor(nday - (nyr-1)*d1yr), '%3.3i');
    hour_num = num2str(mod(nday - d1yr, 1)*24, '%2.2i');
    
    %------------------------------- layers
    for ik = 1:nk
        % layer 
        klay = layers(ik);
        fprintf('\nklay: %s ...\n', num2str(klay));
        %
        %-- uvdp
        dir_path = dir_al{1};
        read_uvdp_GSH;
        fprintf('\nRead UVDP from: %s...\n', file_name);
        fland = dpm <= 1.e-12;
            
        for icel = 1:ncel % [1 2 3 4]
            
            %------------ read uvdp & c
            
            %-- c
            NTRACR = carryTracer_al{icel};
            E = E_al{icel};
            diag_tracer;
            if trac_fid > -1
                tracer(fland) = NaN;
                fprintf(1,'Read tracer-%d from: %s\n', NTRACR, file_name);
            else
                tracer = NaN*zeros(JDM,IDM);
                fprintf(1,'Tracer does NOT exist: %s\n', file_name);
            end
%             c_al{icel}(:,:,ik,it) = tracer;

            %------------ COV matrix
            [varxx, varxy, varyy, xcent, ycent] = tracer_dispersion(tracer,xmesh,ymesh); % tracer.*scp2
            % eigenvalues of the COV matrix
            [phi,lmd1,lmd2] = coordrot_symcomp(varxx,varxy,varxy,varyy);
            %------------ assign
            intc_al{icel}(ik,it) = int_c;
            varxx_al{icel}(ik,it) = varxx;
            varxy_al{icel}(ik,it) = varxy;
            varyy_al{icel}(ik,it) = varyy;
            lmd1_al{icel}(ik,it) = lmd1;
            lmd2_al{icel}(ik,it) = lmd2;
            phi_al{icel}(ik,it) = phi;
        end  % icel   
    end  % ik
end

% save('a_cdispersion.mat','-regexp','c_al');  % COMPASS_dispersion_c

% save('a_cdispersion.mat','-regexp','\w*_al\>'); 

%% plot tracer cloud snapshots, overlapped with initial cloud
% 2-by-2
% load('c_disper_snap.mat')

%--- for initial patch ploted in circle
[r,x0,y0] = deal(75, 900, 400); % see "tracer_per_cell.m"
% len of semi-axis
Scnum = 47;
[ra, rb] = deal(1/Scnum);

% ------
clim = [1e-2 1e0];

cmname = '-haline';
titlstr = {'$\textrm{MEAN}$', '$\textrm{EXP-}\textbf{K}_{red}$',...
    '$\textrm{EXP-}\kappa \mbox{\boldmath $\chi$} $', '$\textrm{FULL}$'};

% titlstr = {'$\textrm{MEAN}$', '$\textrm{EXP-}K_{iso}$',...
%     '$\textrm{EXP-}\kappa \mbox{\boldmath $\chi$} $', '$K_{iso},\,\nabla{K_{iso}}$', '$\textrm{FULL}$'};

% titlstr = {'$\textrm{MEAN}$', '$\textrm{EXP-}\textbf{K}_{iso}$',...
%     '$\textrm{EXP-}\kappa \vec{\chi}$', '$\textrm{FULL}$'};
ids = [1 3 4 6];
it = 1;
ik = 1;
plt_fields = c_zsm3_al(ids); % c_al c_zsm3_al
ncel = numel(plt_fields);

% -------
font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
% [ha, ~] = tight_subplot(2,2,[.02 .02],[.02 .06],[.06 .1]);
[ha, ~] = tight_subplot(2,2,[.02 .02],[.02 .06],[.06 .1]);

for icel = 1:ncel
    axes(ha(icel))
%     subplot(2,2,icel)
    f_do = plt_fields{icel}(:,:,ik,it);
    f_do(f_do < 1e-3) = NaN;
%     f_do = smooth_geom_HYCOM(f_do,scp2,101,101);
    plot_field_model(f_do,plon1d,plat1d,cmname);
    % center of ellipses (K's grid)
    [llx,lly] = m_ll2xy(plon(y0,x0),plat(y0,x0),'clip','off');
    cmap = cmocean(cmname);
    colormap(cmap);
    %----- overlapped with initial patch
    hold on
    % plot ellipse
    h = ellipse(ra,rb,0,llx,lly,'m'); 
    h.LineWidth = 1;

    title([ '$' titlstr{icel} '$'],'interpreter','latex','FontSize',16);
    if icel == 3
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,...
            'ytick',30:5:45,'linewidth',1.2);
    else
        m_grid('linestyle','none','tickdir','out','xtick',280:5:310,'xticklabels',[],...
            'ytick',30:5:45,'yticklabels',[],'linewidth',1.2);
    end
    caxis(clim);
    ax = gca;
    ax.Title.Interpreter = 'Latex';
    ax.Title.String = titlstr{icel};
    ax.ColorScale = 'log';
end
cb = colorbar;
set(cb,'orientation','vertical','Position',[0.92 0.2 0.02 0.63]);
cb.TickDirection = 'out';
cb.LineWidth = 1; % thickness
% cb.TickLength = 0.015;
% cb.XTick = cbTicks;
% cb.XTickLabel = cbTickLabels;
cbarrow;

% set(gcf,'PaperPositionMode','auto'); print(gcf,'fig_cdisp','-dpng','-r600');
% printpdf(gcf,'fig12','-r600')

%% plot total integral of tracer

legends = {'EXP3KL', 'EXP2', 'MEAN', 'FULL'};
linewidths = {2, 2, 2, 2};
linetyles = {'r-','b-','g-','k-'};

ik = 1:30;
font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);

[ha, ~] = tight_subplot(2,1,[.05 .07],[.08 .05],[.07 .05]);
axes(ha(1));

% major 
for icel = 1:ncel
    f_do = sum(intc_al{icel}(ik,:),1);
    plot(1:nt_al,f_do, linetyles{icel},'LineWidth',linewidths{icel});
    hold on;
%     set(gca,'YLim',[-2 2])
end
legend(legends,'Location','best');
title(['Total integral of tracer Z' num2str(ik,'%02d')])
% time-axis
its = 1:4:30;
set(gca,'XTick',its,'XTicklabel','')
% xlabel('Time [d]')
ylabel('Tracer')

%------- plot tracer cloud dispersion
ik = 18;

legends = {'Major (EXP3KL)', 'Major (EXP2)', 'Major (MEAN)', 'Major (FULL)',...
    'Minor (EXP3KL)', 'Minor (EXP2)', 'Minor (MEAN)', 'Minor (FULL)'};
linewidths = {2, 2, 2, 2};

linetyles = {'r-','b-','g-','k-'};

% figure
axes(ha(2));
% major 
for icel = 1:ncel
    f_do = sum(lmd1_al{icel}(ik,:),1);
    plot(1:nt_al,f_do, linetyles{icel},'LineWidth',linewidths{icel});
    hold on;
%     set(gca,'YLim',[-2 2])
end
% minor
linetyles = {'r--','b--','g--','k--'};
for icel = 1:ncel
    f_do = sum(lmd2_al{icel}(ik,:),1);
    plot(1:nt_al,f_do, linetyles{icel},'LineWidth',linewidths{icel});
    hold on;
end
lgd = legend(legends,'Interpreter','latex','FontSize',8);
lgd.Location = 'best'; lgd.Box = 'on';
title(['Tracer dispersion Z' num2str(ik,'%02d')])
% time-axis
relative_times = t_al - t_al(1);
its = 1:4:30;
set(gca,'XTick',its,'XTicklabel',relative_times(its))
xlabel('Time [d]')
ylabel('Dispersion [m^2]')

%% plot dispersion curves for PAPER 

icel_do = [1 2 4 6];
% [lmd1_plt, lmd2_plt] = deal(lmdT1_al, lmdN1_al);
% [lmd1_plt, lmd2_plt] = deal(lmdT_alz, lmdN_alz);
[lmd1_plt, lmd2_plt] = deal(lmdT3_al, lmdN3_al);
% '$\textrm{EXP-}\textbf{K}_{red}$'
% [lmd1_plt, lmd2_plt] = deal(lmdT1_al, lmdN1_al);
if_pltfit = 0;

ik = 1;
firstDay = 22;
relative_times = t_al - firstDay ; % 
its = 1:5:nt_al;

legends = {'Major (MEAN)', 'Major ($\textrm{EXP-}\textbf{K}_{red}$)', 'Major ($\textrm{EXP-}\textbf{K}_{red}$)',...
    'Major ($\textrm{EXP-}\kappa \mbox{\boldmath $\chi$}$)', 'Major (EXP-ADV)', 'Major (FULL)'};
% legends = {'Major (MEAN)', 'Major (FULL)'};
linewidths = {2, 2, 2, 2, 2, 2};
linetyles = {'k--','b-','r-','r-', 'g-','k-'};

% 
font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);

[ha, ~] = tight_subplot(2,1,[.05 .07],[.10 .10],[.10 .10]);
% major 
axes(ha(1))
for icel = icel_do
    f_do = sum(lmd1_plt{icel}(ik,:),1); % loglog(f_do, linetyles{icel},'LineWidth',linewidths{icel});
    if icel == 4

    end
    plot(1:nt_al,f_do, linetyles{icel},'LineWidth',linewidths{icel});
    hold on
end
% 
ax = gca;
ax.TickLength = [.01, .01];
ax.LineWidth = 1.0;
ax.XLim = [1 nt_al+1];
ax.XTick = its;
ax.XTickLabel = '';
% ax.XAxis.Label.Interpreter = 'Latex';
% ax.XAxis.Label.String = 'Day';
ax.YAxis.Label.String = 'Dispersion [m^{2}]';
lgd = legend(legends(icel_do),'Interpreter','latex','FontSize',8);
lgd.Location = 'northwest'; lgd.Box = 'on';
% 
% 
%------ fitted line
if if_pltfit
    dt_in_s = (t_al(2) - t_al(1)) * 24 * 3600; % t_al ~[d]
    linetyles_fit = {'g-','g-','g-'};
    x_text = {19, 22, 15}; y_text = {4e10, 11e10, 11e10};
    incre = {-5e9, 5e9, 5e9};
    for icel = icel_do
        f_do = sum(lmd1_plt{icel}(ik,:),1); % loglog(f_do, linetyles{icel},'LineWidth',linewidths{icel});
        % fit
        it_fit = 15:24; % relative_times(it_fit);
        p = polyfit(it_fit,f_do(it_fit),1); round( p(1)/dt_in_s ,-1)
        f_fitted = polyval(p,it_fit) + incre{icel};
        % plot
        plot(it_fit,f_fitted, linetyles_fit{icel},'LineWidth',1);
        text(x_text{icel},y_text{icel},['k=' num2str(round( p(1)/dt_in_s ,-2),'%4d') 'm^2/s'],'FontSize',10)
        hold on
    end
    hold off
end
%-------------------- minor 
% legends = {'Minor (MEAN)', 'Minor (EXP-KL)', 'Minor (FULL)'};
% legends = {'Minor (MEAN)', 'Minor (FULL)'};

legends = {'Minor (MEAN)', 'Minor ($\textrm{EXP-}\textbf{K}_{red}$)', 'Minor ($\textrm{EXP-}\textbf{K}_{red}$)',...
    'Minor ($\textrm{EXP-}\kappa \mbox{\boldmath $\chi$}$)', 'Minor (EXP-ADV)', 'Minor (FULL)'};
axes(ha(2))
for icel = icel_do
    f_do = sum(lmd2_plt{icel}(ik,:),1); %loglog(f_do, linetyles{icel},'LineWidth',linewidths{icel});
    if icel == 4
    end
    plot(1:nt_al,f_do, linetyles{icel},'LineWidth',linewidths{icel});
    hold on
end
% 
ax = gca;
ax.TickLength = [.01, .01];
ax.LineWidth = 1.0;
ax.XLim = [1 nt_al+1];
% ax.YLim = [0 2e11];
ax.XTick = its;
ax.XTickLabel = relative_times(its);
ax.XAxis.Label.String = 'Day';
ax.YAxis.Label.String = 'Dispersion [m^{2}]';
% 
lgd = legend(legends(icel_do),'Interpreter','latex','FontSize',8);
lgd.Location = 'northwest'; lgd.Box = 'on';
%------ fitted line
if if_pltfit
    dt_in_s = (t_al(2) - t_al(1)) * 24 * 3600; % t_al ~[d]
    linetyles_fit = {'g-','g-','g-'};
    x_text = {21, 23, 18}; y_text = {1.7e10, 2.5e10, 3.3e10};
    incre = {-1e9, 1e9, 1e9};
    for icel = icel_do
        f_do = sum(lmd2_plt{icel}(ik,:),1); % loglog(f_do, linetyles{icel},'LineWidth',linewidths{icel});
        % fit
        it_fit = 18:28; % relative_times(it_fit);
        p = polyfit(it_fit,f_do(it_fit),1); round( p(1)/dt_in_s ,-1)
        f_fitted = polyval(p,it_fit) + incre{icel};
        % plot
        plot(it_fit,f_fitted, linetyles_fit{icel},'LineWidth',1);
        text(x_text{icel},y_text{icel},['k=' num2str(round( p(1)/dt_in_s ,-2),'%4d') 'm^2/s'],'FontSize',10)
        hold on
    end
    hold off
end
% set(gcf,'PaperPositionMode','auto'); print(gcf,'fig14','-dpng','-r600')
% printpdf(gcf,'fig13','-r600')

%% plot intial patch

fland = isnan(c_al{3}(:,:,1));

figure;
[ha, ~] = tight_subplot(1,3,[.04 .01],[.03 .05],[.05 .07]);

c_init = tracer_per_cell(JDM,IDM,12);
c_init(fland) = NaN;
axes(ha(3))
plot_field_model(c_init,plon1d,plat1d,'thermal')
title('Initial tracer patch','interpreter','latex','fontsize',12);

caxis([0 3])
cb = colorbar;
set(cb,'Location','EastOutside','Position',[0.945 0.36 0.015 0.264])

% set(gcf,'PaperPositionMode','auto');
% print(gcf,'COMPASS_disp_cinit','-dpng','-r500')
