% 
% Calc the 'spreading' of F_div (or K) among different tracer pairs, i.e., 
% quantify the dependence of F_div (or K) on tracers.
% 
% bsub -J std -P ome -o pt.o%J -e pt.e%J -W 0:20 -q general -n 1 -R "rusage[mem=25000]" matlab -r fluctuation_K_allterms
% 

addpath(genpath('/nethome/yxl1496/HYCOM'));
addpath(genpath('/nethome/yxl1496/mytoolbox'));

clear
hycom_domain = 'GSH';
read_HYCOM_grid;

%% params for reading

%------------------------------------- input by hand
klay = 24;
wichSM = 4;
ifExtra = 0;
DivAdv = 2; % 1 div; 2 adv
ndist = 3;
ifiso = 1;

%----
[day_s, day_e, dt_save] = deal(100, 100, 2); % 531
t_al = day_s:dt_save:day_e;
nt_al = length(t_al);
d1yr = 365; 

%------------------------------------------
%-------- all tracer sets
carry_al = [1 2 3 4 5 6 7 8 9 10];
trac_comb = nchoosek(carry_al,ndist);
ncomb = size(trac_comb,1);

ncomb_choose = 10;
icomb_choose = sort( randperm(ncomb, ncomb_choose) );

%-------- which filter scales
smdeg_al = [.25 .5 .75 1.0 1.25 1.5]; % half of win leng
dx_mod = 0.02;           
win_len = round(2 * smdeg_al(wichSM) / dx_mod) + 1;

%-------- if smooth the LHS
if ifExtra == 1
    extraStr = 'extraSM';
else
    extraStr = 'noextraSM';
end

%-------- div or adv on the lhs
if DivAdv == 1
    lhsStr = 'div';
elseif DivAdv == 2
    lhsStr = 'adv';
end

%-------- iso model or aniso model
if ifiso
    isoStr = 'iso';
else
    isoStr = 'ani';
end

%------------------------------------- dir
root_dir = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/trac_analysis';
% K dir
termStrs = {'uecs_uece', 'uecs', 'uece', 'usce'};
% termStrs = {'alleddy_cht', 'uecs_cht', 'usce_cht', 'uece_cht'};
% termStrs = {'uscs_dcshs'};
% termStrs = {'usc_dchs'};

ncel = numel(termStrs);
KLdir_al = cell(1,ncel);
for icel = 1:ncel
    KLdir_al{icel} = [root_dir '/diffus_lambda_' lhsStr '_' isoStr '/sm'...
    num2str(win_len,'%03d') '/' termStrs{icel} '/' extraStr ...
    '/Z' num2str(klay,'%02d')];
    fprintf(1,'KL will be read from: %s\n',KLdir_al{icel});
end

%%
flucK_al = cell(1,ncel);
flucK_al(:) = {zeros(JDM,IDM,nt_al)};
flucL_al = cell(1,ncel);
flucL_al(:) = {zeros(JDM,IDM,nt_al)};
% [xx_spm,xy_spm,yx_spm,yy_spm] = deal(zeros(nt,npairs));

%----------------------------------------- loop over time to save memory
for it = 1:nt_al%
    
    %---------------------------------------- current time 
    nday = t_al(it);
    nyr = floor((nday - 1)/d1yr) + 1;   % start from year '0001
    yrStr = num2str(nyr, '%4.4i');   % used for tracer
    dyStr = num2str(floor(nday - (nyr-1)*d1yr), '%3.3i');
    hrStr = num2str(mod(nday - d1yr, 1)*24, '%2.2i');
    disp(['nday =  ' num2str(nday)])
            
    %------------------------------ read K from diff tracer sets
    if ifiso == 1
        fprintf('\n  KL-iso ...\n');
        [K11_al,lmdu_al,lmdv_al] = deal(cell(1,ncel));
        [K11_al(:),lmdu_al(:),lmdv_al(:)] = deal({zeros(JDM,IDM,ncomb_choose)});
    else
        fprintf('\n  KL-aniso ...\n');
        [K11_al,K12_al,K22_al,lmdu_al,lmdv_al] = deal(cell(1,ncel));
        [K11_al(:),K12_al(:),K22_al(:),lmdu_al(:),lmdv_al(:)] = deal({zeros(JDM,IDM,ncomb_choose)});
    end
         
    
    %----------------------- tracer pairs
    for icomb = 1:ncomb_choose
        
        carries = trac_comb(icomb_choose(icomb),:);
        fprintf('\nTracer set: %s...\n',mat2str(carries));

        %-------------------- read 
        for icel = 1:ncel
            fnm = [KLdir_al{icel} '/C' num2str(carries,'%02d') '/K_C',...
                num2str(carries,'%02d') '_D' dyStr 'H' hrStr '_' hycom_domain '.mat'];
            disp(['Reading KL from:  ' fnm])
            struc_rd = load(fnm);
            if ifiso == 1
                [K11,lmdu,lmdv] = deal(struc_rd.K11,struc_rd.lmdu,struc_rd.lmdv);
                K11 = filter_extreme(K11,1,99);
                lmdu = filter_extreme(lmdu,1,99);
                lmdv = filter_extreme(lmdv,1,99);
                %
                [K11_al{icel}(:,:,icomb),lmdu_al{icel}(:,:,icomb),...
                    lmdv_al{icel}(:,:,icomb)] = deal(K11,lmdu,lmdv);
            else
                [K11,K12,K22,lmdu,lmdv] = deal(struc_rd.K11,struc_rd.K12,...
                    struc_rd.K22,struc_rd.lmdu,struc_rd.lmdv);
            end
        end
        
        % eigenvalues of S
%         [phi0_al,kt0_al,kn0_al] = eigens_S(Kxx,Kxy,Kyx,Kyy);
    end
    
    %------------------------------ calc spreading at current time
    dim = 3;
    for icel = 1:ncel
        flucK_al{icel}(:,:,it) = abs(std(K11_al{icel}, 1, dim, 'omitnan') ./ nanmean(K11_al{icel}, dim));
        flucL_al{icel}(:,:,it) = ( abs(std(lmdu_al{icel}, 1, dim, 'omitnan') ./ nanmean(lmdu_al{icel}, dim)) + ...
            abs(std(lmdv_al{icel}, 1, dim, 'omitnan') ./ nanmean(lmdv_al{icel}, dim)) ) ./ 2;
    end
end

savename = ['new_fluc_KL' lhsStr '_' extraStr '_' isoStr '_Z' num2str(klay,'%02d') '.mat'];
save(savename, 'flucK_al','flucL_al','KLdir_al','icomb_choose'); % ,'\<fluc_\w*','\w*0\>'

% save(['fluc_KL' lhsStr '_' extraStr '_' isoStr '.mat'],'flucK_al','flucL_al'); % ,'\<fluc_\w*','\w*0\>'


%% plot snapshots of spreading of Ks from different terms

%
it = 1;
cmapstr = 'amp';
clim = [0 10];

plt_fields = flucL_al;
ncel = numel(plt_fields);

figure
[ha, ~] = tight_subplot(2,2,[.04 .01],[.03 .05],[.05 .07]);
for icel = 1:ncel
    axes(ha(icel));
    f_do = plt_fields{icel}(:,:,it);
    
    %         h = pcolor(plon, plat, f_do); set(h,'EdgeColor', 'none');
    %         set(gca,'xlim',plon([1 end]),'ylim',plat([1 end]));
    %
    plot_field_model(f_do,plon1d,plat1d,cmapstr)
    %         hold on; plot(plon1d,plat1d(269)*ones(1,IDM),'r'); hold on; plot(plon1d,plat1d(805)*ones(1,IDM),'r');
    %
    title(['$ ' termStrs{icel} ... %' , \,Z' num2str(layers(ik),'%02d')
        ' \,\,  [',num2str(nanmean(f_do(:)),'%5.4f'), ';\,\,',...
        num2str(min(f_do(:)),'%5.4f'), ', \,' num2str(max(f_do(:)),'%5.4f')...
        ']$'], 'interpreter','latex','fontsize',12);
    cmocean(cmapstr)
    caxis(clim)
end

set(ha(1:ncel),'XTickLabel',[],'YTickLabel',[]);
cb = colorbar;
%     set(cb,'Location','SouthOutside','Position',[0.30 0.30 0.40 0.03])
set(cb,'Location','EastOutside','Position',[0.945 0.36 0.015 0.264])

%% plot error bars of the spreading
%
it = 1;
plt_fields = flucL_al;
ncel = numel(plt_fields);

figure
[ha, ~] = tight_subplot(2,2,[.05 .05],[.05 .05],[.05 .05]);
for icel = 1:ncel
    axes(ha(icel));
    f_do = plt_fields{icel}(:,:,it);

    [f_m,f_med,f_h,f_l] = stats_ptc(f_do,75,25);

    % prepare data for plotting the error bars
    x = 1 : nt_al;
    y = f_med(x);
    [neg, pos] = deal(y - f_l(x), f_h(x) - y);

    % plot
    axes(ha(icel));
    h = errorbar(x,y,neg,pos);
    ax = gca;
    
    % set properties
    h.Marker = '.';
    h.MarkerSize = 16;
    h.LineStyle = 'none';
    %
    ax.YLim = [0 15];
%     ax.YTick = y(1:2:end);
%     ax.YTickLabel = y(1:2:end);
    ax.Title.Interpreter = 'latex';
    ax.Title.FontSize = 14;

end

%% plot DIST of the spreading
%
it = 1;
plt_fields = flucK_al;
plt_fields2 = flucL_al;
legends = {'K', 'Lambda'};
ncel = numel(plt_fields);

figure
[ha, ~] = tight_subplot(2,2,[.10 .10],[.05 .05],[.10 .05]);
xlim = [0 30];
dx = 1e-1;
edges = xlim(1):dx:xlim(2);

for icel = 1:ncel
    axes(ha(icel));
%     1d, no nan
    f_do = plt_fields{icel}(:,:,it); 
    f_do(isnan(f_do)) = [];
    [histValues,~] = histcounts(f_do,edges,'Normalization','pdf');
    plot(edges(1:end-1),histValues,'-','Color','k','LineWidth',1.2)
    hold on 
    %
    f_do = plt_fields2{icel}(:,:,it);
    f_do(isnan(f_do)) = [];
    [histValues,~] = histcounts(f_do,edges,'Normalization','pdf');
    plot(edges(1:end-1),histValues,'--','Color','k','LineWidth',1.2)
    
    ax = gca;
    ax.TickLength = [.015, .015];
    ax.LineWidth = 1.5;
    ax.YLim = [0 .20];
    ax.Title.String = 'd<c><h>/dt + div<U><c>'; % termStrs{icel}
    ax.Title.String = '<h>d<c>/dt + <U>del<c>'; % termStrs{icel}
%     ax.Title.Interpreter = 'latex';
%     ax.XTick = xlim(1):dx_plt:xlim(2);
    legend(legends,'Location','northeast','FontSize',8);

end

%%
%-------  error bars in ONE plot !
legends = {'alleddy', '$ \textbf{u}^\prime \cdot\nabla \langle{c}\rangle $',...
    '$\langle{\textbf{u}}\rangle \cdot\nabla c^\prime$', ...
    '$\textbf{u}^\prime \cdot\nabla c^\prime$'};
it = 1;
plt_fields = flucK_al;
ncel = numel(plt_fields);
% prepare plot-data
for icel = 1:ncel
    f_do = plt_fields{icel}(:,:,it);
    [f_m(icel),f_med(icel),f_h(icel),f_l(icel)] = stats_ptc(f_do,75,25);
end
% plot
figure
x = 1 : ncel;
y = f_med(x);
[neg, pos] = deal(y - f_l(x), f_h(x) - y);
h = errorbar(x-.1,y,neg,pos); 
ax = gca;
% set properties
h.Marker = '.';
h.Color = 'k';
h.MarkerSize = 16;
h.LineStyle = 'none';
ax.XLim = [0.5 ncel+.5];
ax.XTick = 1:ncel;
ax.XTickLabel = legends;
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0 20];
%     ax.YTick = y(1:2:end);
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

% ---------- for the other param
% prepare plot-data
plt_fields2 = flucL_al;
for icel = 1:ncel
    f_do = plt_fields2{icel}(:,:,it);
    [f_m(icel),f_med(icel),f_h(icel),f_l(icel)] = stats_ptc(f_do,75,25);
end
% plot
x = 1 : ncel;
y = f_med(x);
[neg, pos] = deal(y - f_l(x), f_h(x) - y);
hold on
h = errorbar(x+.1,y,neg,pos); 
ax = gca;
% set properties
h.Marker = '.';
h.Color = 'k';
h.MarkerSize = 16;
h.LineStyle = 'none';
h.Bar.LineStyle = 'dashed';
ax.XLim = [0.5 ncel+.5];
ax.YLim = [0 20];
%     ax.YTick = y(1:2:end);
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

% set(gcf,'PaperPositionMode','auto'); print(gcf,'a','-dpng','-r400');

%% plot norm of K tensor vs. time
[plt_alleddy,plt_uecs,plt_usce,plt_uece] = deal(K0_norm_m,K1_norm_m,K2_norm_m,K3_norm_m);
% [plt_alleddy,plt_uecs,plt_usce,plt_uece] = deal(K0_norm_med,K1_norm_med,K2_norm_med,K3_norm_med);
% [plt_alleddy,plt_uecs,plt_usce,plt_uece] = deal(K0_norm_l,K1_norm_l,K2_norm_l,K3_norm_l);

figure
subplot(211)
plot(plt_alleddy,'.-k','LineWidth',1.2); hold on;
plot(plt_uecs,'.-b','LineWidth',1.2); 
title('domain-averaged norm of K tensor')
set(gca,'XLim',[1 nt])
legend('alleddy','u''<c>','Location','best','fontsize',12);
subplot(212)
plot(plt_usce,'.-g','LineWidth',1.2); hold on;
plot(plt_uece,'.-m','LineWidth',1.2);
set(gca,'XLim',[1 nt])
legend('<u>c''','u''c''','Location','best','fontsize',12);

% set(gcf,'PaperPositionMode','auto'); print(gcf,'Nak_all','-dpng','-r400');


%%
% [plt11,plt12,plt21,plt22] = deal(fluc_K0,fluc_K1,fluc_K2,fluc_K3);
% ylim = [0 10; 0 10];
    
% [plt11,plt12,plt21,plt22] = deal(fluc_lmd0,fluc_lmd1,fluc_lmd2,fluc_lmd3);
% [plt11,plt12,plt21,plt22] = deal(fluc_kt0(:,:,it),fluc_kt1(:,:,it),fluc_kt2(:,:,it),fluc_kt3(:,:,it));
[plt11,plt12,plt21,plt22] = deal(fluc_phi0(:,:,it),fluc_phi1(:,:,it),fluc_phi2(:,:,it),fluc_phi3(:,:,it));

ylim = [0 2; 0 2];

dt_plt = 1;
nt_plt = nsm;
tt_plt = 1:dt_plt:nt_plt;

figure

data2d = reshape(abs(plt11(:,:,tt_plt)),[JDM*IDM nt_plt]);
subplot(221)
boxplot(data2d,'PlotStyle','compact','Whisker',0,'OutlierSize',0.01);
ax = gca;
ax.XTick = 1:dt_plt:nt_plt;
ax.XTickLabel = 2*smdeg_al;
ax.XLim = [0 nt_plt+1];
ax.YLim = ylim(1,:);
title('spreading of Lambda (alleddy)')

data2d = reshape(abs(plt12(:,:,tt_plt)),[JDM*IDM nt_plt]);
subplot(222)
boxplot(data2d,'PlotStyle','compact','Whisker',0,'OutlierSize',0.01);
ax = gca;
ax.XTick = 1:dt_plt:nt_plt;
ax.XTickLabel = 2*smdeg_al;
ax.XLim = [0 nt_plt+1];
ax.YLim = ylim(1,:);
title('u''<c>')

data2d = reshape(abs(plt21(:,:,tt_plt)),[JDM*IDM nt_plt]);
subplot(223)
boxplot(data2d,'PlotStyle','compact','Whisker',0,'OutlierSize',0.01);
ax = gca;
ax.XTick = 1:dt_plt:nt_plt;
ax.XTickLabel = 2*smdeg_al;
ax.XLim = [0 nt_plt+1];
ax.YLim = ylim(2,:);
title('<u>c''')

data2d = reshape(abs(plt22(:,:,tt_plt)),[JDM*IDM nt_plt]);
subplot(224)
boxplot(data2d,'PlotStyle','compact','Whisker',0,'OutlierSize',0.01);
ax = gca;
ax.XTick = 1:dt_plt:nt_plt;
ax.XTickLabel = 2*smdeg_al;
ax.XLim = [0 nt_plt+1];
ax.YLim = ylim(2,:);
title('u''c''')

%% plot K from all tracer sets

cmname = 'balance';
clim = [-1e5 1e5];
plt_fields = Kxx_al{1};

figure
[ha, ~] = tight_subplot(2,2,[.04 .02],[.12 .04],[.02 .05]);

for iplt = 1:1
    axes(ha(iplt));
    
    icomb = iplt + 8;
    f_do = plt_fields(:,:,icomb);
    carries = trac_comb(icomb,:);
    h = pcolor(plon, plat, f_do); set(h,'EdgeColor', 'none');
    set(gca,'xlim',plon([1 end]),'ylim',plat([1 end]));
    title(['$ C (',num2str(carries,'%02d'),'), \,  \lambda_m = ',...
        num2str(nanmean(f_do,'all'),'%5.0f'),...
        '\,\ |\lambda|_m = $',num2str(nanmean(abs(f_do),'all'),'%5.0f')],...
        'interpreter','latex','fontsize',8);
    cmocean(cmname)
    caxis(clim)
end

set(ha(1:4),'XTickLabel',[],'YTickLabel',[]);
cb = colorbar;
set(cb,'Location','SouthOutside','Position',[0.28 0.09 0.40 0.03])

%%
figure
for it = 1:2%:nt
    
    % data
    data = data2d(:,it);
    data(isnan(data)) = [];
    
    % create a different axes for each group:
    ax = axes;
    boxplot(ax,data,'PlotStyle','compact','Whisker',2,'OutlierSize',0.01);
%     ax.XTick = 1.5:2:(groups*2-0.5);
%     ax.XTickLabel = {'a','b','c','v','f'};
    ax.YLim = [-1 1e1];
    box off
    if it == 1 
        ylabel('std/mean') % only for the most right axes 
    else
        ax.YTick = [];
    end
    ax.Position = [corner(it) 0.11 width 0.8];
    
end
axis off

%% plot both K-tensor and K&L !!!

struc_Ktensor = load('fluc_K_div_noextraSM.mat');
struc_KL = load('fluc_KLadv_noextraSM_iso.mat');

flucKtens_al = struc_Ktensor.flucK_al;
flucK_al = struc_KL.flucK_al;
flucL_al = struc_KL.flucL_al;

%-------  error bars in ONE plot !
legends_Kt = {'alleddy', '$\textbf{u}^\prime \langle{c}\rangle$',...
    '$\langle{\textbf{u}}\rangle c^\prime$', '$\textbf{u}^\prime  c^\prime$'};
legends_KL = {'alleddy', '$ \textbf{u}^\prime \cdot\nabla \langle{c}\rangle $',...
    '$\langle{\textbf{u}}\rangle \cdot\nabla c^\prime$', ...
    '$\textbf{u}^\prime \cdot\nabla c^\prime$'};

it = 1;
% plot
figure
[ha, ~] = tight_subplot(1,2,[.04 .10],[.07 .07],[.07 .07]);

% --------
plt_fields = flucKtens_al;
ncel = numel(plt_fields);
% prepare plot-data
for icel = 1:ncel
    f_do = plt_fields{icel}(:,:,it);
    [f_m(icel),f_med(icel),f_h(icel),f_l(icel)] = stats_ptc(f_do,75,25);
end
% plot
axes(ha(1));
x = 1 : ncel;
y = f_med(x);
[neg, pos] = deal(y - f_l(x), f_h(x) - y);
h = errorbar(x,y,neg,pos); 
ax = gca;
% set properties
h.Marker = '.';
h.Color = 'k';
h.MarkerSize = 16;
h.LineStyle = 'none';
ax.XLim = [0.5 ncel+.5];
ax.XTick = 1:ncel;
ax.XTickLabel = legends_Kt;
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0 20];
%     ax.YTick = y(1:2:end);
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
% title
ax.Title.FontSize = 14;
ax.Title.Interpreter = 'latex';
ax.Title.String = '$\textbf{K}$';


% -----------
plt_fields = flucK_al;
ncel = numel(plt_fields);
% prepare plot-data
for icel = 1:ncel
    f_do = plt_fields{icel}(:,:,it);
    [f_m(icel),f_med(icel),f_h(icel),f_l(icel)] = stats_ptc(f_do,75,25);
end
axes(ha(2));

x = 1 : ncel;
y = f_med(x);
[neg, pos] = deal(y - f_l(x), f_h(x) - y);
h = errorbar(x-.1,y,neg,pos); 
ax = gca;
% set properties
h.Marker = '.';
h.Color = 'k';
h.MarkerSize = 16;
h.LineStyle = 'none';
ax.XLim = [0.5 ncel+.5];
ax.XTick = 1:ncel;
ax.XTickLabel = legends_KL;
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0 20];
%     ax.YTick = y(1:2:end);
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

% ---------- for the other param
% prepare plot-data
plt_fields = flucL_al;
for icel = 1:ncel
    f_do = plt_fields{icel}(:,:,it);
    [f_m(icel),f_med(icel),f_h(icel),f_l(icel)] = stats_ptc(f_do,75,25);
end
% plot
x = 1 : ncel;
y = f_med(x);
[neg, pos] = deal(y - f_l(x), f_h(x) - y);
hold on
h = errorbar(x+.1,y,neg,pos); 
ax = gca;
% set properties
h.Marker = '.';
h.Color = 'k';
h.MarkerSize = 16;
h.LineStyle = 'none';
h.Bar.LineStyle = 'dashed';
% title
ax.Title.FontSize = 14;
ax.Title.Interpreter = 'latex';
ax.Title.String = 'K \& $\vec{\lambda}$';
% set(gcf,'PaperPositionMode','auto'); print(gcf,'a','-dpng','-r400');

%}

%%
function [phi, Kt, Kn] = eigens_S(Kxx,Kxy,Kyx,Kyy)

[Kxx_c, Kxy_c, Kyx_c, Kyy_c] = deal( NaN * zeros(size(Kxx)) );
Kxx_c(:,1:end-1) = (Kxx(:,2:end,:) + Kxx(:,1:end-1)) / 2;
Kxy_c(:,1:end-1) = (Kxy(:,2:end) + Kxy(:,1:end-1)) / 2;
Kyx_c(1:end-1,:) = (Kyx(2:end,:) + Kyx(1:end-1,:)) / 2;
Kyy_c(1:end-1,:) = (Kyy(2:end,:) + Kyy(1:end-1,:)) / 2;
[phi, Kt, Kn] = coordrot_symcomp(Kxx_c, Kxy_c, Kyx_c, Kyy_c);

end

function [sp_m,sp_med,sp_ptil_h,sp_ptil_l] = stats_ptc(K,prc_h,prc_l)

sp_m = squeeze(nanmean((K),[1 2]));
sp_med = squeeze(prctile((K),50,[1 2]));
sp_ptil_h = squeeze(prctile((K),prc_h,[1 2]));
sp_ptil_l = squeeze(prctile((K),prc_l,[1 2]));

end