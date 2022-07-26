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
ifTmean = 0;
klay = 24;
wichSM = 4;
ifExtra = 0;
ifDiv = 1;
ifisoK = 0;
ndist = 2;

%----
[day_s, day_e, dt_save] = deal(100, 100, 2); % 531
t_al = day_s:dt_save:day_e;
nt_al = length(t_al);
d1yr = 365; 

%------------------------------------------
%-------- all tracer sets
% carry_al = [1 2 3 4 5 6 7 8 9 10];
carry_al = [1 2 3 4 5 6 7 8 9 10];

trac_comb = nchoosek(carry_al,ndist);
ncomb = size(trac_comb,1);

ncomb_choose = 10;
icomb_choose = sort( randperm(ncomb, ncomb_choose) );

%-------- which filter scales
smdeg_al = [.25 .5 .75 1.0 1.25 1.5]; % half of win leng
dx_mod = 0.02;           
win_len = round(2 * smdeg_al(wichSM) / dx_mod) + 1;

% extra smooth to LHS, if needed
win_extra = win_len;

%-------- if smooth the LHS
if ifExtra == 1
    extraStr = 'extraSM';
else
    extraStr = 'noextraSM';
end

%-------- if use div comp of eddy flux
if ifDiv
    divStr = '_div';
else
    divStr = '';
end

%-------- if calc iso K or K-tensor
if ifisoK == 1
    isoKStr = '_iso';
elseif ifisoK == 0
    isoKStr = '';
elseif ifisoK == 2
    isoKStr = '_isoA';
end

%------------------------------------- dir
root_dir = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/trac_analysis';
% K dir
termStrs = {'uecs_uece', 'uecs', 'uece', 'usce'};
% termStrs = {'alleddy'};
ncel = numel(termStrs);
Kdir_al = cell(1,ncel);
for icel = 1:ncel
    if ifTmean
        Kdir_al{icel} = [root_dir '/diffusivity' isoKStr divStr '/sm'...
            num2str(win_len,'%03d')  '_tmean/' termStrs{icel} '/' extraStr ...
            '/Z' num2str(klay,'%02d')];
    else
        Kdir_al{icel} = [root_dir '/diffusivity' isoKStr divStr '/sm'...
            num2str(win_len,'%03d') '/' termStrs{icel} '/' extraStr '/Z' ...
            num2str(klay,'%02d')];
    end
    fprintf(1,'K will be read from: %s\n',Kdir_al{icel});
end

% [dAdx, dAdy] = calc_GxGy_HYCOM(A,ones(size(depth)),depth,scux,scuy,scvx,scvy,0);
% [uc_EIV, vc_EIV] = deal(dAdy, -dAdx);

%%
if ifisoK == 0
    flucK_al = cell(1,ncel); flucK_al(:) = {zeros(JDM,IDM,nt_al)};
    flucLmd12_al = cell(1,ncel); flucLmd12_al(:) = {zeros(JDM,IDM,nt_al)};
    flucPhi_al = cell(1,ncel); flucPhi_al(:) = {zeros(JDM,IDM,nt_al)};
    %
    flucA_al = cell(1,ncel); flucA_al(:) = {zeros(JDM,IDM,nt_al)};
    fluc_uvA_al = cell(1,ncel); fluc_uvA_al(:) = {zeros(JDM,IDM,nt_al)};
elseif ifisoK == 2
    flucK_al = cell(1,ncel); flucK_al(:) = {zeros(JDM,IDM,nt_al)};
    flucA_al = cell(1,ncel); flucA_al(:) = {zeros(JDM,IDM,nt_al)};
    fluc_uvA_al = cell(1,ncel); fluc_uvA_al(:) = {zeros(JDM,IDM,nt_al)};
end

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
    if ifisoK == 1
        fprintf('\n   K-iso (scalar) ...\n');
        [iso_al] = deal(cell(1,ncel));
        [iso_al(:)] = deal({zeros(JDM,IDM,ncomb_choose)});
    elseif ifisoK == 0 
        fprintf('\n  K-tensor ...\n');
        [xx_al,xy_al,yx_al,yy_al] = deal(cell(1,ncel));
        [xx_al(:),xy_al(:),yx_al(:),yy_al(:)] = deal({zeros(JDM,IDM,ncomb_choose)});
        % symmetric
        [phi_al,lmd1_al,lmd2_al] = deal(cell(1,ncel));
        [phi_al(:),lmd1_al(:),lmd2_al(:)] = deal({zeros(JDM,IDM,ncomb_choose)});
        % anti-symm
        [A_al,uc_A_al,vc_A_al] = deal(cell(1,ncel));
        [A_al(:),uc_A_al(:),vc_A_al(:)] = deal({zeros(JDM,IDM,ncomb_choose)});
    elseif ifisoK == 2
        fprintf('\n  Kiso+A ...\n');
        [Kiso_al,A_al] = deal(cell(1,ncel));
        [Kiso_al(:),A_al(:)] = deal({zeros(JDM,IDM,ncomb_choose)});
        %
        [uc_A_al,vc_A_al] = deal(cell(1,ncel));
        [uc_A_al(:),vc_A_al(:)] = deal({zeros(JDM,IDM,ncomb_choose)});
    end
         
    %----------------------- tracer pairs
    for icomb = 1:ncomb_choose
        
        carries = trac_comb(icomb_choose(icomb),:);
        fprintf('\nTracer set: %s...\n',mat2str(carries));

        %-------------------- read 
        for icel = 1:ncel % term
            fnm = [Kdir_al{icel} '/C' num2str(carries,'%02d') '/K_C',...
                num2str(carries,'%02d') '_D' dyStr 'H' hrStr '_' hycom_domain '.mat'];
            disp(['Reading K from:  ' fnm])
            struc_rd = load(fnm);
            if ifisoK == 1
                Kiso = struc_rd.Kiso;
%                 Kiso = filter_extreme(Kiso,1,99);
                iso_al{icel}(:,:,icomb) = Kiso;
                
            elseif ifisoK == 0
                [Kxx,Kxy,Kyx,Kyy] = deal(struc_rd.Kxx,struc_rd.Kxy,struc_rd.Kyx,struc_rd.Kyy);
                % eigens of S
                [phi, lmd1, lmd2] = eigens_S(Kxx,Kxy,Kyx,Kyy);
                % adv vel of A
                [Kxx_c, Kxy_c, Kyx_c, Kyy_c] = deal( NaN * zeros(size(Kxx)) );
                Kxx_c(:,1:end-1) = (Kxx(:,2:end) + Kxx(:,1:end-1)) / 2;
                Kxy_c(:,1:end-1) = (Kxy(:,2:end) + Kxy(:,1:end-1)) / 2;
                Kyx_c(1:end-1,:) = (Kyx(2:end,:) + Kyx(1:end-1,:)) / 2;
                Kyy_c(1:end-1,:) = (Kyy(2:end,:) + Kyy(1:end-1,:)) / 2;
                A = Kxy_c - Kyx_c;
                [dAdx, dAdy] = calc_GxGy_HYCOM(A,ones(size(depth)),depth,scux,scuy,scvx,scvy,0);
                [uc_EIV, vc_EIV] = deal(dAdy, -dAdx);

%                 Kxx = filter_extreme(Kxx,1,99);
%                 Kxy = filter_extreme(Kxy,1,99);
%                 Kyx = filter_extreme(Kyx,1,99);
%                 Kyy = filter_extreme(Kyy,1,99);
                % A_al,uc_A_al,vc_A_al
                [xx_al{icel}(:,:,icomb),xy_al{icel}(:,:,icomb),...
                    yx_al{icel}(:,:,icomb),yy_al{icel}(:,:,icomb)] = ...
                    deal(Kxx,Kxy,Kyx,Kyy);
                [phi_al{icel}(:,:,icomb),lmd1_al{icel}(:,:,icomb),lmd2_al{icel}(:,:,icomb)] = ...
                    deal(phi, lmd1, lmd2);
                [A_al{icel}(:,:,icomb),uc_A_al{icel}(:,:,icomb),vc_A_al{icel}(:,:,icomb)] = ...
                    deal(A, uc_EIV, vc_EIV);
                
            elseif ifisoK == 2
                [Kiso,A] = deal(struc_rd.Kiso,struc_rd.A);
                %
                [dAdx, dAdy] = calc_GxGy_HYCOM(A,ones(size(depth)),depth,scux,scuy,scvx,scvy,0);
                [uc_EIV, vc_EIV] = deal(dAdy, -dAdx);
                %
                Kiso_al{icel}(:,:,icomb) = Kiso;
                [A_al{icel}(:,:,icomb),uc_A_al{icel}(:,:,icomb),vc_A_al{icel}(:,:,icomb)] = ...
                    deal(A, uc_EIV, vc_EIV);
            end
        end
    end
    
    %------------------------------ calc spreading at current time
    dim = 3;
    % loop over each term!
    for icel = 1:ncel
        if ifisoK == 1
            flucK_al{icel}(:,:,it) = abs( std(iso_al{icel}, 1, dim, 'omitnan') ./ nanmean(iso_al{icel}, dim) );
        elseif ifisoK == 0
            flucK_al{icel}(:,:,it) = ( abs(std(xx_al{icel}, 1, dim, 'omitnan') ./ nanmean(xx_al{icel}, dim)) ...
                + abs(std(xy_al{icel}, 1, dim, 'omitnan') ./ nanmean(xy_al{icel}, dim)) ...
                + abs(std(yx_al{icel}, 1, dim, 'omitnan') ./ nanmean(yx_al{icel}, dim)) ...
                + abs(std(yy_al{icel}, 1, dim, 'omitnan') ./ nanmean(yy_al{icel}, dim)) ) ./ 4;
            % theta, lmd1, lmd2
            flucPhi_al{icel}(:,:,it) = abs(std(phi_al{icel}, 1, dim, 'omitnan') ./ nanmean(phi_al{icel}, dim));
            
            temp12 = (lmd1_al{icel} + lmd2_al{icel}) / 2;
            flucLmd12_al{icel}(:,:,it) = abs(std(temp12, 1, dim, 'omitnan') ./ nanmean(temp12, dim));
            
%             flucLmd12_al{icel}(:,:,it) = ( abs(std(lmd1_al{icel}, 1, dim, 'omitnan') ./ nanmean(lmd1_al{icel}, dim)) + ...
%                 abs(std(lmd2_al{icel}, 1, dim, 'omitnan') ./ nanmean(lmd2_al{icel}, dim)) ) / 2;
            % A, uc*, vc*
            flucA_al{icel}(:,:,it) = abs(std(A_al{icel}, 1, dim, 'omitnan') ./ nanmean(A_al{icel}, dim));
            fluc_uvA_al{icel}(:,:,it) = ( abs(std(uc_A_al{icel}, 1, dim, 'omitnan') ./ nanmean(uc_A_al{icel}, dim)) + ...
                abs(std(vc_A_al{icel}, 1, dim, 'omitnan') ./ nanmean(vc_A_al{icel}, dim)) ) /2;
        elseif ifisoK == 2
            flucK_al{icel}(:,:,it) = abs( std(Kiso_al{icel}, 1, dim, 'omitnan') ./ nanmean(Kiso_al{icel}, dim) );
            % A, uc*, vc*
            flucA_al{icel}(:,:,it) = abs(std(A_al{icel}, 1, dim, 'omitnan') ./ nanmean(A_al{icel}, dim));
            fluc_uvA_al{icel}(:,:,it) = ( abs(std(uc_A_al{icel}, 1, dim, 'omitnan') ./ nanmean(uc_A_al{icel}, dim)) + ...
                abs(std(vc_A_al{icel}, 1, dim, 'omitnan') ./ nanmean(vc_A_al{icel}, dim)) ) /2;
        end
    end
end
%  '_tmean' num2str(ifTmean)
savename = ['new2_fluc_K' isoKStr divStr '_' extraStr '_Z' num2str(klay,'%02d') '.mat'];
save(savename,'-regexp','\<fluc\w*','carry_al','Kdir_al','isoKStr','termStrs',...
    'icomb_choose','trac_comb'); % ,'\<fluc_\w*','\w*0\>'

%
%% plot snapshots of spreading of Ks from different terms

plt_fields = flucL_al;
legends_plt1 = termStrs;

ncel = numel(plt_fields) ; % !!!!!!!!!!!! ---------
% prepare plot-data
for icel = 1:ncel
    % ---- exchange
    f_do = plt_fields{icel}(:,:,it);
    %
    [f_m(icel),f_med(icel),f_h(icel),f_l(icel)] = stats_ptc(f_do,80,10);
end
% plot
figure
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
ax.XTickLabel = legends_plt1; ax.XAxis.TickLabelInterpreter = 'tex';
% ax.TickLabelInterpreter = 'latex';
ax.YLim = ylim;
%     ax.YTick = y(1:2:end);
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
% title
ax.Title.FontSize = 16;
ax.Title.Interpreter = 'latex';
ax.Title.String = '$\textbf{K}$';
% 
ax.XLabel.String = xlabelStr; ax.YLabel.String = ylabelStr;

%{
it = 1;
cmapstr = 'amp';
clim = [0 15];

plt_fields = flucK_al;
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
    %         title(['$ ' titlstr1{icel} ... %' , \,Z' num2str(layers(ik),'%02d')
    %             ' \,\,  [',num2str(nanmean(f_do(:)),'%5.4f'), ';\,\,',...
    %             num2str(min(f_do(:)),'%5.4f'), ', \,' num2str(max(f_do(:)),'%5.4f')...
    %             ']$'], 'interpreter','latex','fontsize',12);
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
plt_fields = flucK_al;
ncel = numel(plt_fields);

figure
[ha, ~] = tight_subplot(2,2,[.05 .05],[.05 .05],[.05 .05]);
for icel = 1:ncel
    axes(ha(icel));
    f_do = plt_fields{icel}(:,:,it);

    [f_m,f_med,f_h,f_l] = stats_ptc(f_do,80,20);

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
    ax.YLim = [0 25];
%     ax.YTick = y(1:2:end);
%     ax.YTickLabel = y(1:2:end);
    ax.Title.Interpreter = 'latex';
    ax.Title.FontSize = 14;

end
%%
%-------  error bars in ONE plot !
legends = {'alleddy', '$u^\prime \langle c \rangle$', '$\langle u \rangle c^\prime$', '$u^\prime c^\prime$'};
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
h = errorbar(x,y,neg,pos); 
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

% set(gcf,'PaperPositionMode','auto'); print(gcf,'a','-dpng','-r400');

%% plot DIST of the spreading
%
it = 1;
plt_fields = flucK_al;
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
%     plot
    [histValues,~] = histcounts(f_do,edges,'Normalization','pdf');
    plot(edges(1:end-1),histValues,'-','Color','k','LineWidth',1.2)
    ax = gca;
    ax.TickLength = [.015, .015];
    ax.LineWidth = 1.5;
    ax.YLim = [0 .20];
    ax.Title.String = termStrs{icel};
%     ax.XTick = xlim(1):dx_plt:xlim(2);
end

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

%}



