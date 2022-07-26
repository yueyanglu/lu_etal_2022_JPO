function [psi,u_psi,v_psi,phi,u_phi,v_phi,output] = uv_decomp(u,v,cxy,alpha,opt)
%Computes the streamfunction and potential and the related rotational and
%  divergent components of a vector field. This function waives the setting
%  of initial guess and the size of variables.
% 

%% check the size of inputs

if ~isequal(length(u),length(v),length(cxy.cux),length(cxy.cvx))
    error('Size of the inputs must be the same !!!');
end

[JDM,IDM] = size(u);
fland = isnan(u);

%% prepare the shape of vars to satisfy 'uv2psiphi.m'

% CHOOSE subregion, with the W/E edges of u and S/N edges of v are deleted
[jjq, iiq] = deal(1:JDM, 1:IDM);                % E.g. 6-6
[jjp, iip] = deal(jjq(1:end-1), iiq(1:end-1));  % E.g. 5-5
[jju, iiu] = deal(jjq(1:end-1), iiq(2:end-1));  % E.g. 5-4
[jjv, iiv] = deal(jjq(2:end-1), iiq(1:end-1));  % E.g. 4-5

% vars
u = u(jju,iiu);
v = v(jjv,iiv);
[cqx,cqy] = deal(cxy.cqx(jjq,iiq), cxy.cqy(jjq,iiq)); 
[cpx,cpy] = deal(cxy.cpx(jjp,iip), cxy.cpy(jjp,iip)); 
[cux,cuy] = deal(cxy.cux(jju,iiu), cxy.cuy(jju,iiu)); 
[cvx,cvy] = deal(cxy.cvx(jjv,iiv), cxy.cvy(jjv,iiv)); 

% size of the grid
[njq, niq] = deal(length(jjq),length(iiq));
[njp, nip] = deal(length(jjp),length(iip));
[nju, niu] = deal(length(jju),length(iiu));
[njv, niv] = deal(length(jjv),length(iiv));

cxy = struct('cpx',cpx,'cpy',cpy,'cqx',cqx,'cqy',cqy,'cux',cux,'cuy',cuy,...
    'cvx',cvx,'cvy',cvy);

% disp 
fprintf('Input u/v size %s\n', mat2str([JDM IDM]))
fprintf('Q-point size %s\n', mat2str([njq niq]))
fprintf('P-point size %s\n', mat2str([njp nip]))
fprintf('U-point size %s\n', mat2str([nju niu]))
fprintf('V-point size %s\n', mat2str([njv niv]))

%% do

%---------------------------------------------- initial guess of psi/phi
% psi0 = zeros(njq,niq); % psi= int_v*dx - int_u*dy, phi = int_u*dx + int_v*dy
% int_vdx = cumsum(v ./ cvx, 1, 'omitnan'); int_vdx(end+1,:) = 0;
% int_udy = cumsum(u ./ cuy, 2, 'omitnan'); int_udy(:,end+1) = 0;
% int_udx = cumsum(u ./ cux, 2, 'omitnan'); int_udx(:,end+1) = 0;
% int_vdy = cumsum(v ./ cvy, 1, 'omitnan'); int_vdy(end+1,:) = 0;
% % 
% psi0(2:end,2:end) = int_vdx - int_udy;
% phi0 = int_udx + int_vdy;
% disp('Random initial guess of psi and phi !!')
psi0 = rand(njq,niq); 
phi0 = rand(njp,nip);

%---------------------------------------------- psi/phi that minimize ja
[psi_ot,phi_ot,output] = uv2psiphi(psi0,phi0,u,v,cxy,alpha,opt);


%% put results onto original grid

[psi,u_psi,v_psi,phi,u_phi,v_phi] = deal(NaN * zeros(JDM,IDM));

%---------------------------------- psi and phi
psi(jjq,iiq) = psi_ot;
phi(jjp,iip) = phi_ot;

%---------------------------------- components of u, v
dpsidy = (psi_ot(2:end,2:end-1) - psi_ot(1:end-1,2:end-1)) .* cuy; % u-
dpsidx = (psi_ot(2:end-1,2:end) - psi_ot(2:end-1,1:end-1)) .* cvx; % v-

dphidy = (phi_ot(2:end,:) - phi_ot(1:end-1,:)) .* cvy; % v-
dphidx = (phi_ot(:,2:end) - phi_ot(:,1:end-1)) .* cux; % u-

% optimal rot and div comp
[u_psi(jju,iiu), v_psi(jjv,iiv)] = deal( - dpsidy, dpsidx);
[u_phi(jju,iiu), v_phi(jjv,iiv)] = deal(   dphidx, dphidy);

% set NaNs in vel
u_psi(fland) = NaN; v_psi(fland) = NaN;
u_phi(fland) = NaN; v_phi(fland) = NaN;


