function [psi,phi] = uv2psiphi(lon,lat,u,v,scux,scvy,ZBC,MB,alpha,fac,period)

if period
    [lon,lat,u] = periodify(lon,lat,u);
    [~,~,v] = periodify(lon,lat,v);
end

flag = isnan(u);

u(flag) = 0;
v(flag) = 0;



