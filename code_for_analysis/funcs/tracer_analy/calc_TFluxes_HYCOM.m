function [utx,vty] = calc_TFluxes_HYCOM(u,v,t,method)
%   Calculates lateral advective tracer fluxes
%   SYNTAX: [utx,vty]=calc_TFluxes_HYCOM(u,v,t,method)
%           utx, vtx -----  [m/s*C]
%           u,v, --- velocities
%           t --- tracer
%           method --- '2nd_order' (default) or '4th_order'

if nargin==3
    method = '2nd_order';
    %disp('Default 2nd order advection is used!')
end

[nj,ni]=size(t);
[utx, vty] = deal(NaN * zeros(nj,ni));
    
if strcmpi(method,'2nd_order')
%
%   Calculate 2nd order tracer fluxes
%   
    iic=2:ni;jjc=2:nj;
    utx(jjc,iic)=.5*u(jjc,iic).*(t(jjc,iic)+t(jjc,iic-1));
    vty(jjc,iic)=.5*v(jjc,iic).*(t(jjc,iic)+t(jjc-1,iic));

elseif strcmpi(method,'4th_order')
%
%   Calculate 4th order tracer fluxes
%
    iic=3:ni-1;jjc=1:nj;
    utx(jjc,iic)=.5*u(jjc,iic).*(1.125*t(jjc,iic)+1.125*t(jjc,iic-1)-0.125*t(jjc,iic+1)-0.125*t(jjc,iic-2));
    for j=1:length(jjc)
%        iip=find(u(jjc(j),iic)>=0);utx(j,iip)=utx(j,iip)-u(jjc(j),iic(iip)).*t(jjc(j),iic(iip)-1);
%        iin=find(u(jjc(j),iic)< 0);utx(j,iin)=utx(j,iin)-u(jjc(j),iic(iin)).*t(jjc(j),iic(iin)  );
    end
    iic=1:ni;jjc=3:nj-1;
    vty(jjc,iic)=.5*v(jjc,iic).*(1.125*t(jjc,iic)+1.125*t(jjc-1,iic)-0.125*t(jjc+1,iic)-0.125*t(jjc-2,iic));
    for i=1:length(iic)
%        jjp=find(v(jjc,iic(i))>=0);vty(jjp,i)=vty(jjp,i)-v(jjc(jjp),iic(i)).*t(jjc(jjp)-1,iic(i));
%        jjn=find(v(jjc,iic(i))< 0);vty(jjn,i)=vty(jjn,i)-v(jjc(jjn),iic(i)).*t(jjc(jjn)  ,iic(i));
    end
end

