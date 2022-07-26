function [K11,K12,K22,lmdu,lmdv] = flx2KLambda_divadv(lhs,del2c,cxx,cxy,cyy,cxp,cyp,ifiso)
% 
% 'flx2KLambda_flxdiv.m' calculates the instantaneous eddy diffusivity 
%  and 'advective vector' from the divergence/advection of eddy tracer flux.
% 
%   Input:
%      lhs   [nj,ni,ndist]: [c*m/s]      
%      del2c [nj,ni,ndist]: [c/m]
%      cxx   [nj,ni,ndist]: 2nd derivative on p-point [c/m]
%      cyy   [nj,ni,ndist]: [c/m]
%      cxp   [nj,ni,ndist]: x-grad of c on p-point [c]
%      cyp   [nj,ni,ndist]: y-grad of c on p-point [c]
%      ifiso              : if use isotropic diffusion (0 for anisotropic)
%
%   Output:
%      K    [nj,ni]: [m2/s]
%      lmd  [nj,ni]: on p-point [m/s]
% 

%% check size of input

% must be 3D
if numel(size(lhs)) >  3
    error('f2K function accepts 3D fields only (j-i-ndist) !!');
end

% inputs must have the same size
if ~isequal(numel(lhs),numel(cxx),numel(cyy),numel(cxp),numel(cyp))
    error('Size of the inputs must be the same !!!');
end

% size
[nj,ni,ndist] = size(lhs);

% iso diffusion: K*laplacian, 3 unknowns, >= 3 tracers;
% anis diffusion: s11*cxx + s12*cxy + s22*cyy, 5 unknowns, >= 5 tracers;
if ifiso
    ndist_min = 3;
else
    ndist_min = 5;
end

if ndist < ndist_min
    error(['Must have at least ' num2str(ndist_min) ' traceres because'...
        ' ifiso = ' num2str(ifiso) ' !!!']);
end

%% calc K

[K11, K12, K22, lmdu, lmdv] = deal(NaN * zeros(nj,ni)); 

% Note this loop will miss the N and E boundaries of the NE
% corner bin.
% ~1min for one snapshot (sequel)
for j = 1:nj
    for i = 1:ni
        
        %------------------------------------ LHS, column vec [div1;div2;div3;..]
        D = squeeze(lhs(j,i,:)); 
        
        %------------------------------------ matrix of Laplacian & grad
        % div_del_c col vec
        Lapl = squeeze(del2c(j,i,:));
        % cxx, cxy, cyy
        Gxx = squeeze(cxx(j,i,:));
        Gxy = squeeze(cxy(j,i,:));
        Gyy = squeeze(cyy(j,i,:));
        % dCdx & dCdy on p
        Gxp = squeeze(cxp(j,i,:));
        Gyp = squeeze(cyp(j,i,:));
        
        % Form the LHS matrix
        % if iso, ndist-by-3; if aniso, ndist-by-5
        if ifiso
            G = [Lapl, -Gxp, -Gyp];
        else
            G = [Gxx, 2*Gxy, Gyy, -Gxp, -Gyp];
        end
        
        %------------------------------------ inverse, -G * KL = D
        if all(~isnan(G),'all') && all(~isnan(D),'all')
            
            if rank(G) == ndist_min
                KL = - G \ D;  % K = - G \ F; or K = - pinv(G) * F
            else
                KL(1:ndist_min) = NaN;
            end
            
            if ifiso
                K11(j,i) = KL(1); K12(j,i) = KL(1); K22(j,i) = KL(1);
                lmdu(j,i) = KL(2); lmdv(j,i) = KL(3);
            else
                K11(j,i) = KL(1); K12(j,i) = KL(2); K22(j,i) = KL(3);
                lmdu(j,i) = KL(4); lmdv(j,i) = KL(5);
            end
        end
    end
end
