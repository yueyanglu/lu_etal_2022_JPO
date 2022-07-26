
% 
%   bsub -J div -P ome -o pt.o%J -e pt.e%J -W 3:00 -q general -n 1 -R "rusage[mem=18000]" matlab -r grad_kappa_chi
% 
addpath(genpath('/nethome/yxl1496/HYCOM'));
addpath(genpath('/nethome/yxl1496/mytoolbox'));

clear
hycom_domain = 'GSH';
read_HYCOM_grid

%% read kappa and chi

[plow, phigh] = deal(1, 99);
klay = 24;
carries = [1 3 5 7 9];
wichTerm = 13;
t_al = 71:1:80; %  71:1:80  71 for dpm
% 
nt_al = length(t_al);
d1yr = 365;
ifiso = 0; % 1 iso; 0 aniso for KL
% 
if wichTerm == 0
    termStr = 'alleddy';
elseif wichTerm == 1
    termStr = 'uecs';
elseif wichTerm == 2
    termStr = 'usce';
elseif wichTerm == 3
    termStr = 'uece';
elseif wichTerm == 13
    termStr = 'uecs_uece';
end
% iso model or aniso model
if ifiso
    isoStr = 'iso';
else
    isoStr = 'ani';
end

% forcing flds dir
forc_dir = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/UVDP_sm101';
KL_dir = ['/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/trac_analysis' ...
    '/diffus_lambda_adv_' isoStr '/sm101/' termStr '/noextraSM/Z' num2str(klay,'%02d') '/C' ...
    num2str(carries,'%02d')];

%------
[K11_tm,K12_tm,K22_tm,lmdu_tm,lmdv_tm] = deal(zeros(JDM,IDM));
dpm_tm = zeros(JDM,IDM);

for it = 1:nt_al 
    
    nday = t_al(it);
    nyr = floor((nday - 1)/d1yr) + 1;   % start from year '0001
    yrStr = num2str(nyr, '%4.4i');   % used for tracer
    dyStr = num2str(floor(nday - (nyr-1)*d1yr), '%3.3i');
    hrStr = num2str(mod(nday - d1yr, 1)*24, '%2.2i');
    fprintf('\nnday: %s ...\n', num2str(nday));

    KL_fnm = [KL_dir '/K_C',num2str(carries,'%02d'),...
        '_D',dyStr,'H',hrStr,'_',hycom_domain,'.mat']; 
    
    % read layer thickness
    [~,~,dpm,~,~] = read_uvdp_GSH_func(forc_dir,dyStr,hrStr,klay,0,0,1,0);
    fprintf('  Read UVDP.\n');
    dpm_tm = dpm_tm + dpm;
    
    % k & Lambda
    struc_load = load(KL_fnm); disp(KL_fnm);
    K11 = struc_load.K11; lmdu = struc_load.lmdu; lmdv = struc_load.lmdv;
    if ifiso
        K12 = NaN * K11;  K22 = NaN * K11;
    else
        K12 = struc_load.K12; K22 = struc_load.K22;
    end
    %filter
    K11 = filter_extreme(K11,plow,phigh);
    K12 = filter_extreme(K12,plow,phigh);
    K22 = filter_extreme(K22,plow,phigh);
    lmdu = filter_extreme(lmdu,plow,phigh);
    lmdv = filter_extreme(lmdv,plow,phigh);

    %
    [K11_tm,K12_tm,K22_tm,lmdu_tm,lmdv_tm] = deal(K11_tm + K11,...
        K12_tm + K12, K22_tm + K22, lmdu_tm + lmdu, lmdv_tm + lmdv); 

end

[K11_tm,K12_tm,K22_tm,lmdu_tm,lmdv_tm] = deal(K11_tm/nt_al,K12_tm/nt_al,K22_tm/nt_al,lmdu_tm/nt_al,lmdv_tm/nt_al);
dpm_tm = dpm_tm/nt_al;


%% 
% By definition, Chi = - del_kappa - del_A + U^xi
% So,  Chi + del_kappa = U^xi - del_A
% 

% kappa & chi*h [m2/s]: p-grid
h = dpm;
[kappa, chiu_p, chiv_p] = deal(K11_tm, lmdu_tm.*h, lmdv_tm.*h);

% del_<h>*kappa [m2/s]: u/v- grid
cflg = 0; 
[delu,delv] = calc_GxGy_HYCOM(kappa,h,depth,scux,scuy,scvx,scvy,cflg);

% chi*h onto u/v- grid
[chiu,chiv] = p2uv(chiu_p,chiv_p);

% Chi + del_kappa (= U^xi - del_A)  [m2/s] : u/v- grid
[resu, resv] = deal(chiu + delu, chiv + delv);

% % curl [F/len = F/m] and div [F*len/area = F/m]
% [fu, fv] = deal(resu, resv); % chiu,chiv   delu_p,delv_p
% curl_chi = calc_curl_HYCOM(fu,fv,scqx,scqy,1);
% div_chi = calc_div_HYCOM(fu,fv,scvx,scuy,scpx.*scpy,1);
% 
% save('chi_delk.mat','kappa','h','chiu_p','chiv_p','resu', 'resv','t_al','klay');

%% calc rot and div comp

%------------------------- params for optmz

% Tikhonov's regularization parameter
alpha = 1e-15;

% options for 'minFunc'
opt.Method = 'lbfgs';
opt.GradObj = 'on';
opt.Display = 'final';
opt.MaxFunEvals = 1e8;
opt.MaxIter = 1e8;
opt.optTol = 1e-8; % 1e-5 (default) Tolerance on the first-order optimality
opt.progTol = 1e-9; % 1e-9 (default) Tolerance on progress in terms of function/parameter changes
opt.DerivativeCheck = 'off';
opt.numDiff = 0;

% mesh size of grid
[cqx,cqy] = deal(1 ./ scqx, 1 ./ scqy); 
[cpx,cpy] = deal(1 ./ scpx, 1 ./ scpy); 
[cux,cuy] = deal(1 ./ scux, 1 ./ scuy); 
[cvx,cvy] = deal(1 ./ scvx, 1 ./ scvy); 
cxy = struct('cpx',cpx,'cpy',cpy,'cqx',cqx,'cqy',cqy,'cux',cux,'cuy',cuy,...
    'cvx',cvx,'cvy',cvy);

%  do
[fu_do,fv_do] = deal(resu, resv);
fprintf('\nCalc div comp at T%f ...\n', nday);
tstart = tic();
[~,uflxR,vflxR,~,uflxD,vflxD,output] = uv_decomp(fu_do,fv_do,cxy,alpha,opt);
cputime = toc(tstart);
toctime = sprintf('Time used is %d seconds', cputime); disp(toctime)

%---- stats

% curl/div of orig flx
curl = calc_curl_HYCOM(fu_do,fv_do,scqx,scqy,1);
div = calc_div_HYCOM(fu_do,fv_do,scvx,scuy,scpx.*scpy,1);

% curl/div of div comp
curl_phi = calc_curl_HYCOM(uflxD,vflxD,scqx,scqy,1);
div_phi = calc_div_HYCOM(uflxD,vflxD,scvx,scuy,scpx.*scpy,1);

% curl/div of rot comp
curl_psi = calc_curl_HYCOM(uflxR,vflxR,scqx,scqy,1);
div_psi = calc_div_HYCOM(uflxR,vflxR,scvx,scuy,scpx.*scpy,1);

% disp
RD_curl = abs( (curl_psi - curl) ./ curl);
RD_div = abs( (div_phi - div) ./ div);

RD_curl_perc = nansum(RD_curl(:) < 0.01) / nansum(~isnan(RD_curl(:)));
RD_div_perc = nansum(RD_div(:) < 0.01) / nansum(~isnan(RD_div(:)));

fprintf('Max and mean of abs(curl_phi) %f, %f \n',...
    max(abs(curl_phi),[],'all'),nanmean(abs(curl_phi),'all'));
fprintf('Max and mean of abs(div_psi) %f, %f \n',...
    max(abs(div_psi),[],'all'),nanmean(abs(div_psi),'all'));
fprintf('Perct of RD btw curl and curl_psi less than .01 is %f \n',RD_curl_perc);
fprintf('Perct of RD btw div and div_phi less than .01 is %f \n\n',RD_div_perc);

% save
save('chi_delk_divrot.mat','uflxR','vflxR','uflxD','vflxD', 'output','opt');


%% plot hist 
%

plt_fields = {lhs1, lhs2};
titls = {'curl(\chi + \langle{h}\rangle\nabla\kappa)', 'div(\chi + \langle{h}\rangle\nabla\kappa)'};

figure
subplot(121)
histogram(plt_fields{1},[-6e-4:5e-5:6e-4],'Normalization','probability')
title([ '$' titls{1} '$'],'interpreter','latex','FontSize',16);
subplot(122)
histogram(plt_fields{2},[-6e-4:5e-5:6e-4],'Normalization','probability')
title([ '$' titls{2} '$'],'interpreter','latex','FontSize',16);

%%

fplt_al = {uc_EIV, vc_EIV};
clim = [-.2, .2];
figure
for icel = 1:2
    subplot(1,2,icel)
    f_do = fplt_al{icel};
    f_do = filter_extreme(f_do,1,99);
    f_do = smooth_geom_HYCOM(f_do,scpx.*scpy,win_extra,win_extra);

    plot_field_model( f_do,plon1d,plat1d,'balance')
    caxis(clim)
    colorbar
end

%% plot
cmname = 'balance';

titls = {'\chi_{u}', '\chi_{v}', '\chi_{u} + \partial_x \kappa', '\chi_{v} + \partial_y \kappa'};
titls = {'\chi_{u}', '\chi_{v}', '-\partial_x \kappa', '-\partial_y \kappa', ...
    '\chi_{u} + \partial_x \kappa', '\chi_{v} + \partial_y \kappa'};
% titls = {'curl(\chi + \langle{h}\rangle\nabla\kappa)', 'div(\chi + \langle{h}\rangle\nabla\kappa)', '\chi_{u} + \partial_x \kappa', '\chi_{v} + \partial_y \kappa'};

% plt_fields = {chi_u, chi_v, delu_p, delv_p, res_u, res_v};  % delu_p, delv_p
plt_fields = {chiu_p, chiv_p, resu, resv};  clim = [-3e2 3e2]; % delu_p, delv_p
% plt_fields = {curl_chi, div_chi};  clim = [-6e-3 6e-3];

% plt_fields = {K11_tm, K12_tm};  clim = [-6e5 6e5];
% plt_fields = {lmdu_tm, lmdv_tm};  clim = [-1e1 1e1];


font = 'DejaVu Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
[ha, ~] = tight_subplot(2,2,[.07 .02],[.05 .08],[.05 .05]);
for icel = 1:4
    axes(ha(icel))
    plot_field_model( plt_fields{icel},plon1d,plat1d,cmname)
    title([ '$' titls{icel} '$'],'interpreter','latex','FontSize',16);
    m_grid('linestyle','none','tickdir','out','xtick',280:5:310,'xticklabels',[],...
        'ytick',30:5:45,'yticklabels',[],'linewidth',1.2);
    caxis(clim)
    %
    cb = colorbar;
    cb.Orientation = 'vertical';
    cb.Title.String = 'm {s^{-1}}';
%     cb.Title.Position = [80 -30 0];
end
%}
