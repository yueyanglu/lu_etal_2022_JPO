clear

%% relative err  btw orig & recons flux div (read data: 'err_reconsflx_fromK.m')

% legends_Kt = {'$\nabla\cdot \left( \textbf{U}^\prime \langle{c}\rangle + \langle{\textbf{U}}\rangle c^\prime + \textbf{U}^\prime  c^\prime \right)$', ...
%     '$\nabla\cdot \textbf{U}^\prime \langle{c}\rangle$',...
%     '$\nabla\cdot \langle{\textbf{U}}\rangle c^\prime$',...
%     '$\nabla\cdot \textbf{U}^\prime  c^\prime$'};
% % 
% legends_KL = {'$ \textbf{U}^\prime \cdot\nabla \langle{c}\rangle +\langle{\textbf{U}}\rangle \cdot\nabla c^\prime + \textbf{U}^\prime \cdot\nabla c^\prime$',...
%     '$ \textbf{U}^\prime \cdot\nabla \langle{c}\rangle $',...
%     '$\langle{\textbf{U}}\rangle \cdot\nabla c^\prime$', ...
%     '$\textbf{U}^\prime \cdot\nabla c^\prime$'};

legends_Kt = {'$\nabla\cdot \left( \textbf{U}^\prime \langle{c}\rangle+ \textbf{U}^\prime  c^\prime \right)$', ...
    '$\nabla\cdot \textbf{U}^\prime \langle{c}\rangle$',...
    '$\nabla\cdot \textbf{U}^\prime  c^\prime$'};

legends_KL = {'$ \textbf{U}^\prime \cdot\nabla \langle{c}\rangle + \textbf{U}^\prime \cdot\nabla c^\prime$',...
%     '$ \textbf{U}^\prime \cdot\nabla \langle{c}\rangle $',...
    '$\textbf{U}^\prime \cdot\nabla c^\prime$'};
legends_str = legends_Kt;

ylim = [0 5];
yticks = ylim(1):1:ylim(2);
% font = 'DejaVu Sans'; figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
figure
[ha, ~] = tight_subplot(3,1,[.07 .05],[.08 .06],[.06 .02]);


% ------------ read data for all four terms
ncel = 3; % # of terms
terms = [13 1 3];
for icel = 1:ncel
    
    %--- read err data
%     fnm = ['errsTm13_K_div_2080_methd2_8samp.mat']; % adv_iso
%     fnm = ['errsTm' num2str(terms(icel),'%01d') '_adv_iso_2080_methd2_8samp.mat'];
    fnm = ['errsTm' num2str(terms(icel),'%01d') '_K_div_2080_methd2_8samp.mat'];
%     fnm = ['errsTm' num2str(terms(icel),'%01d') '_K_2080_methd2.mat'];
    struc_rd = load(fnm);
    
    err_med = struc_rd.err_med; err_m = struc_rd.err_m;
    err_l = struc_rd.err_l; err_h = struc_rd.err_h;
    ncK_al = struc_rd.ncK_al; ncF_al = struc_rd.ncF_al; 
    
    % ---- eleminate NaNs in err data
    id_comb = find(~isnan(err_med));
%     id_comb = icomb_choose;
    err_med_do = err_med(id_comb);
    err_m_do = err_m(id_comb);
    err_l_do = err_l(id_comb);
    err_h_do = err_h(id_comb);
    %
    ncomb_do = length(id_comb);
    ncK_do = ncK_al(id_comb);
    ncF_do = ncF_al(id_comb);

    % ----- data for plot
    dx = 1;
    x = 1:dx:ncomb_do;
    y = err_med_do(x); [neg, pos] = deal(y - err_l_do(x), err_h_do(x) - y);
    xlabels = ncK_do(x);
    
    % ----- plot
    axes(ha(icel));
    h = errorbar(x,y,neg,pos);
    ax = gca;
    % set properties of plot
    h.Marker = '.';
    h.Color = 'k';
    h.MarkerSize = 12;
    h.LineStyle = 'none';
    ax.TickDir = 'in'; 
    ax.XLim = [x(1)-2*dx x(end)+2*dx];
    ax.YLim = ylim;
    ax.XTick = x;
    ax.XTickLabel = xlabels; % ncF_al
    ax.YTick = yticks;

    ax.Title.Interpreter = 'Latex';
    ax.Title.String = legends_str{icel}; % titleStr
    ax.Title.FontSize = 16;
    ax.TickLength = [.005, .005];
    ax.LineWidth = 1.0;
    %
    ax.TickLabelInterpreter = 'Latex'; 
    ax.XAxis.FontSize = 10; ax.XAxis.TickLength = [.01 .01];
    ax.YAxis.FontSize = 12; ax.YAxis.TickLength = [.01 .01];
    ax.XLabel.Interpreter = 'Latex'; ax.XAxis.Label.String = '$N$'; ax.XAxis.Label.FontSize = 12;
    ax.YLabel.Interpreter = 'Latex'; ax.YAxis.Label.String = '$\varepsilon$'; ax.YAxis.Label.FontSize = 12;
    
    if icel < ncel
        ax.XTickLabel = '';
        ax.XAxis.Label.String = ''; 
    end
    if icel < ncel
%         ax.YAxis.Label.String = '';
    end
    
    
end

% set(gcf,'PaperPositionMode','auto'); 
% print(gcf,'error_recons_Ktens','-dpng','-r600'); % error_recons_KL error_recons_Ktens
