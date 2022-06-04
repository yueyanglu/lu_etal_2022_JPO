clear
hycom_domain = 'GSH';
read_HYCOM_grid
load('c_dispersion_evolution.mat')

%% plot dispersion curves 

icel_do = [1 2 4 6];
[lmd1_plt, lmd2_plt] = deal(lmdT1_al, lmdN1_al);
% [lmd1_plt, lmd2_plt] = deal(lmdT3_al, lmdN3_al);
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