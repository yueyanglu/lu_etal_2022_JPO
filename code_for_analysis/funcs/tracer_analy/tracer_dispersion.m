function [varxx, varxy, varyy, xcent, ycent] = tracer_dispersion(c,xmesh,ymesh)
% 
% c (2d): can be tracer or weighted tracer
% xmesh [m]: cumsum(scux,2) - scux(:,1)
% 
% Garrett, 1983; Wagner et al., 2019
% 

%------------ position of center of mass (CoM) at current layer
%-- total int of tracer
int_c = sum(c, 'all', 'omitnan');

%-- total int of tracer weighted by position
int_cx = sum(c .* xmesh , 'all', 'omitnan');
int_cy = sum(c .* ymesh , 'all', 'omitnan');

%-- postion of CoM, same unit with 'xmesh'
[xcent, ycent] = deal( int_cx / int_c, int_cy / int_c );

%------------ second moment (dispersion / variance with respect to CoM)
varxx = sum(c .* (xmesh - xcent) .* (xmesh - xcent), 'all', 'omitnan') / int_c;
varxy = sum(c .* (xmesh - xcent) .* (ymesh - ycent), 'all', 'omitnan') / int_c;
varyy = sum(c .* (ymesh - ycent) .* (ymesh - ycent), 'all', 'omitnan') / int_c;
