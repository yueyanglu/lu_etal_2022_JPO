function [Gx,Gy] = calc_GxGy_HYCOM(c,dpm,depth,scux,scuy,scvx,scvy,flg)
%
% Calculate grid-weighted tracer gradients for the diffusion operator.
%   Gradients [c*m/m] can be multiplied by the grid length, consistent with
%   the HYCOM use of diffusivity. 
% 
% Outputs:
%   Gx,Gy -- gradients on u/v points, same size with c. The unit is 
%            controled by 'flg' 
% 
% Inputs:
%     c -- tracer fld
%   dpm -- layer thickness [m]
%   depth -- depth [m]
%   scux,scuy,scvx,scvy - 
%   flg -- flag for whether multiply gradient by grid len. '0' is NO [c],
%          others are YES [m*c].
% 
% Syntax: [Gx,Gy] = calc_GxGy_HYCOM(c,dpm,depth,scux,scuy,scvx,scvy,flg)
% 

% -------------------------------------------- params
epsil = 1.e-11; 
aspmax = 2.0; 
onemu = 1.e-12; % thickness of vanishing layer
[nj,ni] = size(c);

% -------------------------------------------- check the sizes
ss1 = size(dpm)==size(c);  ss2 = size(depth)==size(c);
ss3 = size(scux)==size(c); ss4 = size(scuy)==size(c);
ss5 = size(scvx)==size(c); ss6 = size(scvy)==size(c);
if prod([ss1 ss2 ss3 ss4 ss5 ss6]) == 0
    error('Inconsistent matrix sizes!');
end

if nargin < 8
    flg = 1;
end

% -------------------------------------------- zonal grad
[jc, ic] = deal(1:nj, 2:ni);

% layer thickness at the interface by harmonic mean [m]
factor = harmon(max(dpm(jc,ic-1),onemu), max(dpm(jc,ic),onemu)); 

% tracer gradient, h*c/dx [c] 
Gx = factor .* (c(jc,ic) - c(jc,ic-1)) ./  max(scux(jc,ic), epsil); 

% grad multiply by grid len [c*m] 
if flg
    % Largest grid spacing (within limits) used in all diffusion coef
    aspux = min(  max(scux(jc,ic), scuy(jc,ic)),...
        min(scux(jc,ic), scuy(jc,ic))*aspmax  );
    Gx = aspux .* Gx;
end

% NaN 
ffl = isnan(depth(jc,ic) .* depth(jc,ic-1));
Gx(ffl) = NaN;

% same size with tracer fld
Gx = [NaN * ones(nj,1), Gx];

% -------------------------------------------- meridional grad
[jc, ic] = deal(2:nj, 1:ni); 

% layer thickness at the interface [m]
factor = harmon( max(dpm(jc-1,ic),onemu), max(dpm(jc,ic),onemu) );

% tracer gradient in [c] 
Gy = factor .* (c(jc,ic) - c(jc-1,ic)) ./  max(scvy(jc,ic), epsil);

% grad multiply by grid len [c*m]  
if flg
    % [m]
    aspvy = min(  max(scvx(jc,ic), scvy(jc,ic)),...
        min(scvx(jc,ic), scvy(jc,ic))*aspmax  );
    Gy = aspvy .* Gy;
end

% NaN 
ffl = isnan(depth(jc,ic) .* depth(jc-1,ic));
Gy(ffl) = NaN;

% size
Gy = [NaN * ones(1,ni); Gy];

end


% harmonic mean of layer thickiness at the interface
function h = harmon(x,y)
h = 2 * x .* y ./ (x + y);
end
