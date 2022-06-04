clear
hycom_domain = 'GSH';
read_HYCOM_grid;
load('stats_tr_evolution.mat')

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