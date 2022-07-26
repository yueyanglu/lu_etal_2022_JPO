
clear
fname = 'csave.mat';
load(fname,'u','v','thck'); 

hycom_domain = 'GSH';
read_HYCOM_grid

% flag of land and zero thck (p-)
fland_p = isnan(depth) | thck <= 1e-12;
[J,I] = find(fland_p == 1);

% flags for u and v, make sure vels surrounding land is NaN. That is, one
% NaN at p-cell yields 2 NaNs in u and 2 in v. This is what HYCOM has done, 
% but we also want treat zero-thickness point as NaNs.
fland_u = fland_p;
fland_v = fland_p;
fland_q = fland_p;

for k = 1:length(I)
    
    % index of u-point at the EAST edge of the land-cell (J(k), I(k))
    [jE,iE] = deal(J(k), I(k)+1);
    % index of v-point at the NORTH edge of land-cell
    [jN,iN] = deal(J(k)+1, I(k));
    % index of q-point at the NORTH-EAST corner of land-cell
    [jNE,iNE] = deal(J(k)+1, I(k)+1);
    
    %
    iE(iE > IDM) = IDM;
    jN(jN > JDM) = JDM;
    jNE(jNE > JDM) = JDM;
    iNE(iNE > IDM) = IDM;

    fland_u(jE,iE) = 1;
    fland_v(jN,iN) = 1;
   
    fland_q(jE,iE) = 1; fland_q(jN,iN) = 1; fland_q(jNE,iNE) = 1;
    
end

% set uv to NaN
u(fland_u) = NaN;
v(fland_v) = NaN;

%% build appropriate grid

% CHOOSE subregion  E.g., 6min for 500-by-500
[jjq, iiq] = deal(1:100, 1:500);          % E.g. 6-6
[jjp, iip] = deal(jjq(1:end-1), iiq(1:end-1));  % E.g. 5-5
[jju, iiu] = deal(jjq(1:end-1), iiq(2:end-1));  % E.g. 5-4
[jjv, iiv] = deal(jjq(2:end-1), iiq(1:end-1));  % E.g. 4-5

%------------------- model vars in the sub-region
thck = thck(jjp,iip);
u = u(jju,iiu);
v = v(jjv,iiv);

fland_q = fland_q(jjq,iiq);
fland_p = fland_p(jjp,iip);
fland_u = fland_u(jju,iiu);
fland_v = fland_v(jjv,iiv);

%------------------- subregion's mesh
[qlon,qlat] = deal(qlon(jjq,iiq),qlat(jjq,iiq));
[plon,plat] = deal(plon(jjp,iip),plat(jjp,iip));
[ulon,ulat] = deal(ulon(jju,iiu),ulat(jju,iiu));       
[vlon,vlat] = deal(vlon(jjv,iiv),vlat(jjv,iiv));  

% size of the grid
[njq, niq] = size(qlon);
[njp, nip] = size(plon);
[nju, niu] = size(ulon);
[njv, niv] = size(vlon);

% mesh size of C-grid
[cqx,cqy] = deal(1 ./ scqx(jjq,iiq), 1 ./ scqy(jjq,iiq)); 
[cpx,cpy] = deal(1 ./ scpx(jjp,iip), 1 ./ scpy(jjp,iip)); 
[cux,cuy] = deal(1 ./ scux(jju,iiu), 1 ./ scuy(jju,iiu)); 
[cvx,cvy] = deal(1 ./ scvx(jjv,iiv), 1 ./ scvy(jjv,iiv)); 

cxy = struct('cpx',cpx,'cpy',cpy,'cqx',cqx,'cqy',cqy,'cux',cux,'cuy',cuy,...
    'cvx',cvx,'cvy',cvy);

%% do

% Tikhionov's regularization parameter
alpha = 1e-16;

% Which optimz function to use. 1 for MATLAB's fminunc; 2 for minFunc
whichmethod = 2; 

%---------------------------------------------- options for optmz

if whichmethod == 1
    opt = optimoptions('fminunc');
    opt.Algorithm = 'quasi-newton'; % trust-region, quasi-newton
    opt.HessUpdate = 'bfgs';
    opt.SpecifyObjectiveGradient = true;
    opt.Display = 'iter';
    opt.MaxFunEvals = 1e6;
    opt.MaxIter = 1e3;
    opt.OptimalityTolerance = 1e-10;
    opt.DerivativeCheck = 'off'; % check if supplied grad match FD approximations
    opt.FiniteDifferenceType = 'central';
else
    opt.Method = 'lbfgs';
    opt.GradObj = 'on';
    opt.Display = 'final';
    opt.MaxFunEvals = 1e6;
    opt.MaxIter = 1e6;
    opt.optTol = 1e-10; % Tolerance on the first-order optimality
    opt.DerivativeCheck = 'off';
    opt.numDiff = 0;
end

%---------------------------------------------- initial guess of psi/phi
psi0 = rand(njq,niq);
phi0 = rand(njp,nip);
% mask land with 0 (though not necessary), more land points than A'*x
psi0(fland_q) = 0;
phi0(fland_p) = 0;

%---------------------------------------------- psi/phi that minimize ja
tic;
[psi_ot,phi_ot] = uv2psiphi(psi0,phi0,u,v,cxy,alpha,opt,whichmethod);
toc;

%---------------------------------------------- psi/phi to uv

dpsidy = (psi_ot(2:end,2:end-1) - psi_ot(1:end-1,2:end-1)) .* cuy; % u-
dpsidx = (psi_ot(2:end-1,2:end) - psi_ot(2:end-1,1:end-1)) .* cvx; % v-

dphidy = (phi_ot(2:end,:) - phi_ot(1:end-1,:)) .* cvy; % v-
dphidx = (phi_ot(:,2:end) - phi_ot(:,1:end-1)) .* cux; % u-

% optimal rot and div comp
[u_psi_re, v_psi_re] = deal( - dpsidy, dpsidx);
[u_phi_re, v_phi_re] = deal(   dphidx, dphidy);
[u_re, v_re] = deal(u_psi_re + u_phi_re, v_psi_re + v_phi_re);

