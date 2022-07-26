function divF = calc_div_HYCOM(Fx,Fy,scvx,scuy,scp2,weight_opt)
%
%   Calculate horizontal flux divergence on a HYCOM grid. 
%   Assume that Fx is given on the U-grid, and Fy -- on the V-grid.
%   Fluxes are usually given in terms of layer transports uch (m2/s*c),
%   to get flux divergence units, one needs to divide by dpm later 
%
%   Syntax: divF=calc_div_HYCOM(Fx,Fy,scvx,scuy,scp2,dpm)
% 
%           divF -- p-points
%           Fx -- flux in the x-direction
%           Fy -- flux in the y-direction
%           scp2 = scpx*scpy
%

[njp,nip] = size(Fx);
[njp2,nip2] = size(scp2);

divF = NaN * zeros(njp,nip);

if nargin < 6
    weight_opt = 1;
end

if njp2 ~= njp || nip2 ~= nip
    error('Grid and flux dimensions do not match!');
end

[jc, ic] = deal(1:njp-1, 1:nip-1);

if weight_opt == 1
    divF(jc,ic) = (Fx(jc,ic+1).*scuy(jc,ic+1) - Fx(jc,ic).*scuy(jc,ic) + ...
        Fy(jc+1,ic).*scvx(jc+1,ic) - Fy(jc,ic).*scvx(jc,ic)) ./ scp2(jc,ic); 
elseif weight_opt == 2
    divF(jc,ic) = (Fx(jc,ic+1) - Fx(jc,ic) + Fy(jc+1,ic)- Fy(jc,ic));
end
