
clear
hycom_domain = 'GSH';
read_HYCOM_grid;

%%
% ---- plot diffusivities (left) and advective fluxes (right)

f_Ktens = load('fluc_K_div_noextraSM_Z24.mat'); 
% new2: (lmd1+lmd2)/2
f_Kred = load('fluc_K_isoA_div_noextraSM_Z24.mat');
f_KL = load('fluc_KLadv_noextraSM_iso_Z24.mat');

% diffus
flucLmd12_al = f_Ktens.flucLmd12_al; % flucPhi_al flucLmd12_al
flucKiso_al = f_Kred.flucK_al;
flucKappa_al = f_KL.flucK_al;

% adv
flucuvA_Ktens_al = f_Ktens.fluc_uvA_al; 
flucuvA_Kred_al = f_Kred.flucA_al;
flucChi_al = f_KL.flucL_al;


%% plot

%------- 
% together
ylim = [0 10];
markers = {'*' 'o' '.'};
colors = {'k' 'k' 'k'};
marksz = 12;
xlabelStr = '';
ylabelStr = 'Ratio';
xlabels = {'$\textbf{U}^\prime \left( \langle{c}\rangle + c^\prime \right)$',...
    '$\textbf{U}^\prime \langle{c}\rangle$', '$\textbf{U}^\prime  c^\prime$'};

% left
title_L = '';%'diffusive part'; \, \&\, \theta
legendStr_L = {'$(\lambda_1+\lambda_2)/2$', '$K_{iso}$', '$\kappa$'};
plt_fields_L1 = flucLmd12_al;
plt_fields_L2 = flucKiso_al;
plt_fields_L3 = flucKappa_al;

% right
title_R = '';%'advective part';
legendStr_R = {'$-\hat{\textbf{z}} \times \nabla A$',... % _{\textbf{K}}
    '$-\hat{\textbf{z}} \times \nabla A_{red}$',...
    ' \boldmath{$\chi$} '};
plt_fields_R1 = flucuvA_Ktens_al;
plt_fields_R2 = flucuvA_Kred_al;
plt_fields_R3 = flucChi_al;


% ---------------------------- 
it = 1;
font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
[ha, ~] = tight_subplot(1,2,[.05 .07],[.10 .10],[.08 .03]);

% ----------------------------  plot left panel (three diffus)
dx_plt = 0.2;
ncel = numel(plt_fields_L1) - 1; % !!!!!!!!!!!! ---------

% ------ plot left-1
axes(ha(1));
for icel = 1:ncel
    % ---- exchange
    f_do = plt_fields_L1{icel}(:,:,it);
    %
    [f_m(icel),f_med(icel),f_h(icel),f_l(icel)] = stats_ptc(f_do,75,25);
end
% plot
x = 1 : ncel;
y = f_med(x);
[neg, pos] = deal(y - f_l(x), f_h(x) - y);
h = errorbar(x-dx_plt,y,neg,pos); 
ax = gca;
% set properties
h.Marker = markers{1}; h.Color = colors{1}; h.MarkerSize = marksz;
h.LineStyle = 'none';
ax.XLim = [0.5 ncel+.5];
ax.XTick = 1:ncel;
ax.XTickLabel = xlabels; ax.XAxis.TickLabelInterpreter = 'latex';
% ax.TickLabelInterpreter = 'latex';
ax.YLim = ylim;
%     ax.YTick = y(1:2:end);
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
% title
ax.Title.FontSize = 16;
ax.Title.String = title_L;
% 
ax.XLabel.String = xlabelStr; ax.YLabel.String = ylabelStr;

hold on

% ------ plot left - 2
for icel = 1:ncel
    % ---- exchange
    f_do = plt_fields_L2{icel}(:,:,it);
    %
    [f_m(icel),f_med(icel),f_h(icel),f_l(icel)] = stats_ptc(f_do,75,25);
end
% plot
x = 1 : ncel;
y = f_med(x);
[neg, pos] = deal(y - f_l(x), f_h(x) - y);
h = errorbar(x,y,neg,pos); 
ax = gca;
% set properties
h.Marker = markers{2}; h.Color = colors{2}; h.MarkerSize = marksz;
h.LineStyle = 'none';

hold on

% ------ plot left - 3
for icel = 1:ncel
    % ---- exchange
    f_do = plt_fields_L3{icel}(:,:,it);
    %
    [f_m(icel),f_med(icel),f_h(icel),f_l(icel)] = stats_ptc(f_do,75,25);
end
% plot
x = 1 : ncel;
y = f_med(x);
[neg, pos] = deal(y - f_l(x), f_h(x) - y);
h = errorbar(x+dx_plt,y,neg,pos); 
ax = gca;
% set properties
h.Marker = markers{3}; h.Color = colors{3}; h.MarkerSize = marksz;
h.LineStyle = 'none';

lgd = legend(legendStr_L,'Interpreter','latex');
lgd.Orientation = 'vertical';
lgd.FontSize = 12;

% ------------- right 
% ------ plot r-1
axes(ha(2));
for icel = 1:ncel
    % ---- exchange
    f_do = plt_fields_R1{icel}(:,:,it);
    %
    [f_m(icel),f_med(icel),f_h(icel),f_l(icel)] = stats_ptc(f_do,75,25);
end
% plot
x = 1 : ncel;
y = f_med(x);
[neg, pos] = deal(y - f_l(x), f_h(x) - y);
h = errorbar(x-dx_plt,y,neg,pos); 
ax = gca;
% set properties
h.Marker = markers{1}; h.Color = colors{1}; h.MarkerSize = marksz;
h.LineStyle = 'none';
ax.XLim = [0.5 ncel+.5];
ax.XTick = 1:ncel;
ax.XTickLabel = xlabels; ax.XAxis.TickLabelInterpreter = 'latex';
% ax.TickLabelInterpreter = 'latex';
ax.YLim = ylim;
%     ax.YTick = y(1:2:end);
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
% title
ax.Title.FontSize = 16;
ax.Title.String = title_R;
% 
ax.XLabel.String = xlabelStr; 

hold on

% ------ plot r - 2
for icel = 1:ncel
    % ---- exchange
    f_do = plt_fields_R2{icel}(:,:,it);
    %
    [f_m(icel),f_med(icel),f_h(icel),f_l(icel)] = stats_ptc(f_do,75,25);
end
% plot
x = 1 : ncel;
y = f_med(x);
[neg, pos] = deal(y - f_l(x), f_h(x) - y);
h = errorbar(x,y,neg,pos); 
ax = gca;
% set properties
h.Marker = markers{2}; h.Color = colors{2}; h.MarkerSize = marksz;
h.LineStyle = 'none';

hold on

% ------ plot r - 3
for icel = 1:ncel
    % ---- exchange
    f_do = plt_fields_R3{icel}(:,:,it);
    %
    [f_m(icel),f_med(icel),f_h(icel),f_l(icel)] = stats_ptc(f_do,80,20);
end
% plot
x = 1 : ncel;
y = f_med(x);
[neg, pos] = deal(y - f_l(x), f_h(x) - y);
h = errorbar(x+dx_plt,y,neg,pos); 
ax = gca;
% set properties
h.Marker = markers{3}; h.Color = colors{3}; h.MarkerSize = marksz;
h.LineStyle = 'none';

lgd = legend(legendStr_R,'Interpreter','latex');
lgd.Orientation = 'vertical';
lgd.FontSize = 12;


% -------- plot (a), (b)...
dim1 = [0.09 0.80 0.07 0.075];
dim2 = [0.57 0.80 0.07 0.075];

annotation('textbox',dim1,'String','a','Margin',1,'FitBoxToText','on','EdgeColor','none',...
    'BackgroundColor','w','fontsize',18,'fontname',font,'FontWeight','bold','HorizontalAlignment','center');
annotation('textbox',dim2,'String','b','Margin',1,'FitBoxToText','on','EdgeColor','none',...
    'BackgroundColor','w','fontsize',18,'fontname',font,'FontWeight','bold','HorizontalAlignment','center');


% set(gcf,'PaperPositionMode','auto'); print(gcf,'fig5','-dpng','-r600');
