function [cs,uscs,vscs,uecs,vecs,usce,vsce,uece,vece] = trac2flx_HYCOM(c,u,v,scp2,scu2,scv2,len,methd,ifinterp)
% 
% 'trac2flux_HYCOM.m' calculates the 2D large-scale tracer and eddy fluxes 
% fields from the given tracer field.
% 
%   Input:
%      c    [nj,ni]: instantaneous tracer field
%      u    [nj,ni]: Zonal flx 
%      v    [nj,ni]: Meridional flx
%      len : window length of spatial running mean
%      methd: method for calculating flux
%      ifinterp : if interpolate tracer (used for particle-based method)
%
%   Output:
%      cs : large-scale tracer
%      ucE : zonal eddy flux, uc - <u><c>, with the unit of [c]*[uflx]
%      ucE_S: spatially smoothed eddy flux, < uc - <u><c> >
% 
%  [cs,ucE,vcE,uecs,vecs,usce,vsce,uece,vece] = 
%         trac2flx_HYCOM(c,fu,fv,scp2,scu2,scv2,len,methd,ifinterp)
% 


%% check size of input

% must be 2D
if numel(size(c)) > 2
    error('c2flx function accepts 2D fields only !!');
end

% c,uflx,vflx must be in the same size
if ~isequal(numel(c),numel(u),numel(v))
    error('Size of the inputs must be the same !!!');
end

% 
if ifinterp
    disp('Interp tracer field...');
end

%%

%---------------------------------------- interp trac fld (or NOT)
if ifinterp
    % 2d instantaneous fld
    X = plon(~isnan(c));
    Y = plat(~isnan(c));
    V = c(~isnan(c));
    F = scatteredInterpolant(X,Y,V,'linear','none');
    
    % renew tracer
    c = F(plon,plat);    
    clearvars X Y V F
end

%----------------------------------------------- mean flds, <c> <u> <v>
cs = smooth_geom_HYCOM(c, scp2, len, len);
us = smooth_geom_HYCOM(u, scu2, len, len);
vs = smooth_geom_HYCOM(v, scv2, len, len);

%----------------------------------------------- eddy flds, c' u' v'
ce = c - cs;
ue = u - us;
ve = v - vs;
    
%----------------------------------------------- fluxes

% <u><c>
[uscs, vscs] = calc_TFluxes_HYCOM(us, vs, cs, methd);

% u'<c>
[uecs, vecs] = calc_TFluxes_HYCOM(ue, ve, cs, methd);

% <u>c'
[usce, vsce] = calc_TFluxes_HYCOM(us, vs, ce, methd);
    
% u'c'
[uece, vece] = calc_TFluxes_HYCOM(ue, ve, ce, methd);

