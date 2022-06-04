clear
hycom_domain = 'GSH';
read_HYCOM_grid
load('c_patch_snap.mat')

%% plot tracer cloud snapshots, overlapped with initial cloud
% 2-by-2

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
