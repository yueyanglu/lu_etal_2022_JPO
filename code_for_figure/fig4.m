
clear
load('stat_divadv_noextraSM.mat')

%% magnitude of components of adv/div (for read, 'comp_mag_divadv.m')

ik = 2;
klay = layers(ik);

% ---- div
% ylim1 = {[0 4e-6], [0 4e-6], [0 4e-7] [0 4e-7]};
% ylim1 = {[0 4e-6], [0 1e-4], [0 4e-6] [0 4e-6]};
% ylim1 = {[0 6e-5], [0 3e-3], [0 6e-5] [0 6e-5]}; % 24 noSM
ylim1 = {[0 2e-5], [0 8e-4], [0 2e-5] [0 5e-5]};

legendStr1 = {'$ \nabla\cdot \langle{\textbf{U}}\rangle  \langle{c}\rangle $',...
    '$\nabla\cdot \textbf{U}^\prime   \langle{c}\rangle $',...
    '$\nabla\cdot\langle{\textbf{U}}\rangle   c^\prime$', ...
    '$\nabla\cdot\textbf{U}^\prime   c^\prime$'};
plt1_struc = struct('m',divuc_m_al, 'med',divuc_med_al, 'h',divuc_h_al, 'l',divuc_l_al);
% plt1_struc = struct('m',cdivu_m_al, 'med',cdivu_med_al, 'h',cdivu_h_al, 'l',cdivu_l_al);

% ---- adv
% ylim2 = {[0 6e-5], [0 6e-5], [0 6e-5] [0 6e-5]};
% ylim2 = {[0 2e-6], [0 2e-6], [0 2e-6] [0 2e-6]};
% ylim2 = {[0 4e-7], [0 4e-7], [0 4e-7] [0 4e-7]};
ylim2 = {[0 8e-6], [0 2e-5], [0 2e-5] [0 3e-5]}; % 24 NOSM
% ylim2 = {[0 1e-5], [0 4e-6], [0 4e-6] [0 4e-6]};

legendStr2 = {'$ \langle{\textbf{U}}\rangle \cdot\nabla \langle{c}\rangle $',...
    '$ \textbf{U}^\prime \cdot\nabla \langle{c}\rangle $',...
    '$\langle{\textbf{U}}\rangle \cdot\nabla c^\prime$', ...
    '$\textbf{U}^\prime \cdot\nabla c^\prime$'};
plt2_struc = struct('m',advuc_m_al, 'med',advuc_med_al, 'h',advuc_h_al, 'l',advuc_l_al);


% legendStr1 = {'$ \langle{c}\rangle \nabla\cdot \langle{\textbf{U}}\rangle $',...
%     '$\langle{c}\rangle \nabla\cdot \textbf{U}^\prime $',...
%     '$ c^\prime \nabla\cdot\langle{\textbf{U}}\rangle $', ...
%     '$ c^\prime \nabla\cdot\textbf{U}^\prime  $'};
% plt1_struc = struct('m',cdivu_m_al, 'med',cdivu_med_al, 'h',cdivu_h_al, 'l',cdivu_l_al);
% ylim1 = {[0 2e-5], [0 8e-4], [0 5e-7] [0 1e-5]};


colorStr = {'k', 'k','k','k'};


font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);

[ha, ~] = tight_subplot(2,4,[.07 .06],[.10 .05],[.09 .02]);

% ----------------- plot div

for isubplt = 1:8
    
    axes(ha(isubplt));
    if isubplt <= 4
        icel = isubplt;
        plt_struc = plt1_struc; % div
        ylim = ylim1;
        legendStr = legendStr1;
    else
        icel = isubplt - 4;
        plt_struc = plt2_struc; % adv
        ylim = ylim2;
        legendStr = legendStr2;
    end

    %
    dx = 1;
    x = 1:dx:nt_al;
    y = plt_struc(icel).med(ik,x);
    [neg, pos] = deal(y - plt_struc(icel).l(ik,x), plt_struc(icel).h(ik,x) - y);
    h = errorbar(x,y,neg,pos); hold on
    % text box for legend
    htext = text(.5,.9,legendStr{icel},'Units','normalized','Interpreter','latex');
    htext.EdgeColor = 'k';
    htext.FontSize = 10;

    ax = gca;
    % set properties of plot
    h.Marker = '.';
    h.Color = colorStr{icel};
    h.MarkerSize = 16;
    h.LineStyle = '-';
    ax.XLim = [0.5 nt_al+.5];
    ax.YLim = ylim{icel};
    ax.XTick = x;
    ax.XTickLabel = t_al(x) - t_al(1);
    ax.Title.FontSize = 14;
    ax.TickLength = [.01, .01];
    ax.LineWidth = 1.0;
    %
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 14;
    if isubplt <= 4; ax.XTickLabel = ''; end
    if isubplt <= 4; ax.XTickLabel = ''; end

    if isubplt > 4; ax.XAxis.Label.String = 'Day'; end
    if isubplt == 1 || isubplt == 5; ax.YAxis.Label.String = 'm/s{\cdot}c'; end
    %
end

% set(gcf,'PaperPositionMode','auto'); print(gcf,'mag_divad','-dpng','-r600');
