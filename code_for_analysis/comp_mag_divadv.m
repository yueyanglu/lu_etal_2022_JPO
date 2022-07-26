

clear

% ---- cluster-dependent home dir
homedir = getenv('HOME');
addpath(genpath([homedir '/HYCOM']));
addpath(genpath([homedir '/mytoolbox']));
addpath(genpath([homedir '/MyFuncs']));

% 
hycom_domain = 'GSH';
read_HYCOM_grid;
scp2 = scpx .* scpy;
scu2 = scux .* scuy;
scv2 = scvx .* scvy;

%%
[day_s, day_e, dt_save] = deal(56, 296, 80); % 56, 296, 80   51-671 (T1D051-Y2D306)
t_al = day_s:dt_save:day_e;
nt_al = length(t_al);
d1yr = 365; 
[UFLX_INDEX, VFLX_INDEX, DPM_INDEX, DPI_INDEX] = deal(1);

win_len = 101;
layers = [15 24];
ifExtra = 0;
if ifExtra == 1
    extraStr = 'extraSM';
else
    extraStr = 'noextraSM';
end
nk = length(layers);
carryTracer = 1;
E = '/glade/scratch/yueyanglu/hra_expt/UVDP/out2';
% /projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/expt/ctest4/exp0
dir_path = '/glade/scratch/yueyanglu/hra_expt/UVDP';

[phigh, plow] = deal(80, 20);

%%
ncel = 4;
[divuc_al, cdivu_al, advuc_al] = deal( cell(1,ncel) );
[divuc_al(:),cdivu_al(:),advuc_al(:)] = deal( {zeros(JDM,IDM)} ); % ,nk,nt_al
% 
[divuc_m_al, divuc_med_al, divuc_h_al, divuc_l_al] = deal( cell(1,ncel) );
[cdivu_m_al, cdivu_med_al, cdivu_h_al, cdivu_l_al] = deal( cell(1,ncel) );
[advuc_m_al, advuc_med_al, advuc_h_al, advuc_l_al] = deal( cell(1,ncel) );
[divuc_m_al(:),divuc_med_al(:),divuc_h_al(:),divuc_l_al(:)] = deal( {zeros(nk,nt_al)} );
[cdivu_m_al(:),cdivu_med_al(:),cdivu_h_al(:),cdivu_l_al(:)] = deal( {zeros(nk,nt_al)} );
[advuc_m_al(:),advuc_med_al(:),advuc_h_al(:),advuc_l_al(:)] = deal( {zeros(nk,nt_al)} );

for it = 1:nt_al
    
    %------------------------------- current time
    nday = t_al(it);
    nyr = floor((nday - 1)/d1yr) + 1;   % start from year '0001
    yrStr = num2str(nyr, '%4.4i');   % used for tracer
    dyStr = num2str(floor(nday - (nyr-1)*d1yr), '%3.3i');
    hrStr = num2str(mod(nday - d1yr, 1)*24, '%2.2i');
    
    for ik = 1:nk
        %------------------------------- read uvdp c
        klay = layers(ik);
        fprintf('\nklay: %s ...\n', num2str(klay));
        
        %---- uvdp
        [uflx,vflx,dpm,~,file_name] = read_uvdp_GSH_func(dir_path,dyStr,hrStr,klay,1,1,1,1);
        fprintf('\nRead UVDP from: %s...\n', file_name);
        % <forcings>, NaNs will replace 0s in forcings.
        dpmS = smooth_geom_HYCOM(dpm, scp2, win_len, win_len);
        uflxS = smooth_geom_HYCOM(uflx, scu2, win_len, win_len);
        vflxS = smooth_geom_HYCOM(vflx, scv2, win_len, win_len);
        [uflxE, vflxE] = deal(uflx - uflxS, vflx - vflxS);
        %
        fland = dpm <= 1.e-12;
        
        %---- c
        NTRACR = carryTracer;
        [tracer,file_name] = diag_tracer_func(E,yrStr,dyStr,hrStr,klay,NTRACR);
        fprintf(1,'\nRead tracer-%d from: %s\n', NTRACR, file_name);
        %
        cS = smooth_geom_HYCOM(tracer, scp2, win_len, win_len);
        cE = tracer - cS;
        
        %------------------------------- different terms of c-flx
        [~,uscs,vscs,uecs,vecs,usce,vsce,uece,vece] = trac2flx_HYCOM(...
            tracer,uflx,vflx,scp2,scu2,scv2,win_len,'4th_order',0);
        
        %------------------------------- adv/div of c-flx terms
        %--- for <U><c>
        divuc_al{1} = calc_div_HYCOM(uscs,vscs,scvx,scuy,scp2,1);
        cdivu_al{1} = cS .* calc_div_HYCOM(uflxS,vflxS,scvx,scuy,scp2,1);
        advuc_al{1} = divuc_al{1} - cdivu_al{1};
        
        %--- for U'<c>
        divuc_al{2} = calc_div_HYCOM(uecs,vecs,scvx,scuy,scp2,1);
        cdivu_al{2} = cS .* calc_div_HYCOM(uflxE,vflxE,scvx,scuy,scp2,1);
        advuc_al{2} = divuc_al{2} - cdivu_al{2};
        
        %--- for <U>c'
        divuc_al{3} = calc_div_HYCOM(usce,vsce,scvx,scuy,scp2,1);
        cdivu_al{3} = cE .* calc_div_HYCOM(uflxS,vflxS,scvx,scuy,scp2,1);
        advuc_al{3} = divuc_al{3} - cdivu_al{3};
        
        %--- for U'c'
        divuc_al{4} = calc_div_HYCOM(uece,vece,scvx,scuy,scp2,1);
        cdivu_al{4} = cE .* calc_div_HYCOM(uflxE,vflxE,scvx,scuy,scp2,1);
        advuc_al{4} = divuc_al{4} - cdivu_al{4};
        
        %--- stats !!!
        for icel = 1:ncel
            % div_uc
            f_do = divuc_al{icel};
            if ifExtra
                f_do = smooth_geom_HYCOM(f_do, scp2, win_len, win_len);
                disp('Extra SM on div_uc terms...')
            end
            f_do = abs(f_do);
            [divuc_m_al{icel}(ik,it),divuc_med_al{icel}(ik,it),...
                divuc_h_al{icel}(ik,it),divuc_l_al{icel}(ik,it)] ...
                = stats_ptc(f_do,phigh,plow);
            % cdiv_u
            f_do = cdivu_al{icel};
            if ifExtra
                f_do = smooth_geom_HYCOM(f_do, scp2, win_len, win_len);
                disp('Extra SM on cdiv_u terms...')
            end
            f_do = abs(f_do);
            [cdivu_m_al{icel}(ik,it),cdivu_med_al{icel}(ik,it),...
                cdivu_h_al{icel}(ik,it),cdivu_l_al{icel}(ik,it)] ...
                = stats_ptc(f_do,phigh,plow);
            % adv
            f_do = advuc_al{icel};
            if ifExtra
                f_do = smooth_geom_HYCOM(f_do, scp2, win_len, win_len);
                disp('Extra SM on adv terms...')
            end
            f_do = abs(f_do);
            [advuc_m_al{icel}(ik,it),advuc_med_al{icel}(ik,it),...
                advuc_h_al{icel}(ik,it),advuc_l_al{icel}(ik,it)] ...
                = stats_ptc(f_do,phigh,plow);
        end
        
    end % ik
end % it

savename = [ 'a_newrigt_stat_divadv_' extraStr '.mat'];
save(savename,'-regexp','layers','t_al','nt_al','phigh','plow',...
    '\w*_m_al\>','\w*_med_al\>','\w*_h_al\>','\w*_l_al\>','ifExtra') 
       
%% plot vertical profile of stats 
%{
xlim = [0 1e-3];
% fld to be analyzed
[f_m,f_med,f_h,f_l] = deal(divuc_m_al{1},divuc_med_al{1},divuc_h_al{1},divuc_l_al{1});
% [f_m,f_med,f_h,f_l] = deal(uadvc_m,uadvc_med,uadvc_h,uadvc_l);
% prepare data for plotting the vertical profile
y = 1 : KDM;
x = f_m(y);
[neg, pos] = deal(x - f_l(y), f_h(y) - x);

%
figure
% plot
% axes(ha(icel));
h = errorbar(x,y,neg,pos,'horizontal');
ax = gca;

% set properties
h.Marker = '.';
h.MarkerSize = 16;
h.LineStyle = 'none';
%
ax.YDir = 'reverse';
ax.XLim = xlim;
ax.YLim = [y(1)-1 y(end)+1];
ax.YTick = y(1:2:end);
ax.YTickLabel = y(1:2:end);
ax.Title.Interpreter = 'latex';
ax.Title.FontSize = 14;

%% norm of flx vs. time 
%
figure
[ha, ~] = tight_subplot(2,2,[.04 .10],[.07 .07],[.07 .07]);

[plt1,plt2,plt3,plt4] = deal(normD0_m,normD1_m,normD2_m,normD3_m);
axes(ha(1));
% plot(plt1,'k-','LineWidth',1.2,'MarkerSize',10); hold on;
plot(plt2,'b-','LineWidth',1.2,'MarkerSize',10); 
% title('ROT COMP removed')
set(gca,'XLim',[1 nt],'XTick',1:20:nt,'XTickLabel',0:20:nt-1,'YLim',[5 9])
legend('$ \textbf{u}^\prime\left<c\right> $',...
    'Location','southeast','fontsize',12,'interpreter','latex');
title('$ \left|\left| \mathbf{F_{div}}  \right|\right|_m $','fontsize',16,'interpreter','latex');
ylabel('$C \cdot \mathrm{m^2} \, \mathrm{s^{-1}}$','interpreter','latex')

axes(ha(3));
plot(plt3,'g-','LineWidth',1.2,'MarkerSize',10); hold on;
plot(plt4,'m-','LineWidth',1.2,'MarkerSize',10);
set(gca,'XLim',[1 nt],'XTick',1:20:nt,'XTickLabel',0:20:nt-1,'YLim',[0 1])
legend('$ \left<\textbf{u}\right>c^\prime $',' $ \textbf{u}^\prime c^\prime $',...
    'Location','northeast','fontsize',12,'interpreter','latex');
xlabel('Day','interpreter','latex')
ylabel('$C \cdot \mathrm{m^2} \, \mathrm{s^{-1}}$','interpreter','latex')
set(ha(1),'XTickLabel','')



[plt1,plt2,plt3,plt4] = deal(divM_m,div1_m,div2_m,div3_m);
axes(ha(2));
% plot(plt1,'k-','LineWidth',1.2,'MarkerSize',10); hold on;
plot(plt2,'b-','LineWidth',1.2,'MarkerSize',10); 
% title('ROT COMP removed')
set(gca,'XLim',[1 nt],'XTick',1:20:nt,'XTickLabel',0:20:nt-1,'YLim',[1e-3 2e-3])
legend('$ \textbf{u}^\prime\left<c\right> $',...
    'Location','southeast','fontsize',12,'interpreter','latex');
title('$ \left|\left| \nabla \cdot \mathbf{F}  \right|\right|_m $','fontsize',16,'interpreter','latex');
ylabel('$C \cdot \mathrm{m} \, \mathrm{s^{-1}}$','interpreter','latex')

axes(ha(4));
plot(plt1,'r-','LineWidth',1.2,'MarkerSize',10); hold on;
plot(plt3,'g-','LineWidth',1.2,'MarkerSize',10); hold on;
plot(plt4,'m-','LineWidth',1.2,'MarkerSize',10);
set(gca,'XLim',[1 nt],'XTick',1:20:nt,'XTickLabel',0:20:nt-1,'YLim',[0 2e-4])
legend('$ \left<\textbf{u}\right>\left<c\right> $', ...
    '$ \left<\textbf{u}\right>c^\prime $',' $ \textbf{u}^\prime c^\prime $',...
    'Location','northeast','fontsize',12,'interpreter','latex');
xlabel('Day','interpreter','latex')
ylabel('$C \cdot \mathrm{m} \, \mathrm{s^{-1}}$','interpreter','latex')
set(ha(2),'XTickLabel','')
%


%% stats of div/adv vs. time 

% ---- div
% ylim = {[0 6e-5], [0 3e-3], [0 6e-5] [0 6e-5]};
% legendStr = {'$ \nabla\cdot \langle{\textbf{U}}\rangle  \langle{c}\rangle $',...
%     '$\nabla\cdot \textbf{U}^\prime   \langle{c}\rangle $',...
%     '$\nabla\cdot\langle{\textbf{U}}\rangle   c^\prime$', ...
%     '$\nabla\cdot\textbf{U}^\prime   c^\prime$'};
% plt1_struc = struct('m',divuc_m_al, 'med',divuc_med_al, 'h',divuc_h_al, 'l',divuc_l_al);

% ---- adv
ylim = {[0 6e-5], [0 6e-5], [0 6e-5] [0 6e-5]};
legendStr = {'$ \langle{\textbf{U}}\rangle \cdot\nabla \langle{c}\rangle $',...
    '$ \textbf{U}^\prime \cdot\nabla \langle{c}\rangle $',...
    '$\langle{\textbf{U}}\rangle \cdot\nabla c^\prime$', ...
    '$\textbf{U}^\prime \cdot\nabla c^\prime$'};
plt1_struc = struct('m',advuc_m_al, 'med',advuc_med_al, 'h',advuc_h_al, 'l',advuc_l_al);

colorStr = {'k', 'k','k','k'};

figure
[ha, ~] = tight_subplot(2,2,[.06 .05],[.10 .05],[.10 .10]);

% ----------------- plot div

for icel = 1:4
    
    axes(ha(icel));

    dx = 2;
    x = 1:dx:nt_al;
    y = plt1_struc(icel).med(x);
    [neg, pos] = deal(y - plt1_struc(icel).l(x), plt1_struc(icel).h(x) - y);
    h = errorbar(x,y,neg,pos); hold on
    % text box for legend
    htext = text(.7,.9,legendStr{icel},'Units','normalized','Interpreter','latex');
    htext.EdgeColor = 'k';
    htext.FontSize = 12;
    
    ax = gca;
    % set properties of plot
    h.Marker = '.';
    h.Color = colorStr{icel};
    h.MarkerSize = 16;
    h.LineStyle = '-';
    ax.XLim = [0 nt_al];
    ax.YLim = ylim{icel};
    ax.XTick = x;
    ax.XTickLabel = t_al(x) - t_al(1);
    if icel == 1 || icel == 2; ax.XTickLabel = ''; end
    ax.Title.Interpreter = 'Latex';
    ax.Title.String = titleStr{icel};
    ax.Title.FontSize = 14;
    ax.TickLength = [.005, .005];
    ax.LineWidth = 1.0;
    %
    ax.XAxis.FontSize = 12;
    ax.XAxis.Label.Interpreter = 'latex';
    ax.YAxis.FontSize = 14;
    ax.YAxis.Label.Interpreter = 'latex';
    if icel == 3 || icel == 4; ax.XAxis.Label.String = '$\mathrm{Day}$'; end
    if icel == 1 || icel == 3; ax.YAxis.Label.String = '$\mathrm{m/s \cdot c}$'; end
    %
end

% set(gcf,'PaperPositionMode','auto'); print(gcf,'a','-dpng','-r400');
%}

