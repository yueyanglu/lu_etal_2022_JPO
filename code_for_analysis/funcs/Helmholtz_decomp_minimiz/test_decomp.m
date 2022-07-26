
% bsub -J read -P ome -o de.o%J -e de.e%J  -W 1:00 -q general -n 1 -R "rusage[mem=17000]" matlab -r test_decomp

addpath(genpath('/nethome/yxl1496/HYCOM'));
addpath(genpath('/nethome/yxl1496/mytoolbox'));

%% read flux, 5h20min for one snapshot
it = 22;

hycom_domain = 'GSH';
read_HYCOM_grid

% read <u'c'>
flxname = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/trac_mod_dia/Z24/trac_flx/flx_C01_SM051_GSH.mat';
load(flxname,'ucE_s','vcE_s');

% read h
varmod_file = '/nethome/yxl1496/HYCOM/vars_mod_Z24_24_T140_240.nc';
thck = ncread(varmod_file,'layers_thcks');

% flux to be decomposed  <u'c'> * h
[uch,vch] = calc_TFluxes_HYCOM(ucE_s(:,:,it),vcE_s(:,:,it),thck(:,:,it),...
    '2nd_order');
uch(abs(uch) <= 1.e-12) = NaN;
vch(abs(vch) <= 1.e-12) = NaN;
clearvars ucE_s vcE_s thck

% mesh size of C-grid
[cqx,cqy] = deal(1 ./ scqx, 1 ./ scqy); 
[cpx,cpy] = deal(1 ./ scpx, 1 ./ scpy); 
[cux,cuy] = deal(1 ./ scux, 1 ./ scuy); 
[cvx,cvy] = deal(1 ./ scvx, 1 ./ scvy); 
cxy = struct('cpx',cpx,'cpy',cpy,'cqx',cqx,'cqy',cqy,'cux',cux,'cuy',cuy,...
    'cvx',cvx,'cvy',cvy);

%% do
%------------------------- Tikhionov's regularization parameter
alpha = 1e-15;

%------------------------- options for 'minFunc'
opt.Method = 'lbfgs';
opt.GradObj = 'on';
opt.Display = 'iter';
opt.MaxFunEvals = 1e6;
opt.MaxIter = 1e6;
opt.optTol = 1e-8; % 1e-5 (default) Tolerance on the first-order optimality
opt.progTol = 1e-9; % 1e-9 (default)
opt.DerivativeCheck = 'off';
opt.numDiff = 0;

%------------------------- all output vars are on orig grid
tic;
[psi,u_psi,v_psi,phi,u_phi,v_phi,output] = uv_decomp(uch,vch,cxy,alpha,opt);
toc;

%% save 

save('cc.mat','psi*','phi*','u*','v*','output','opt');

