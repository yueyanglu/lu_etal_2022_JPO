function curlF = calc_curl_HYCOM(Fx,Fy,scqx,scqy,weight_opt)
% 
% curl = dv/dx - du/dy
% 

[njq,niq] = size(Fx);
[njq2,niq2] = size(scqx);
if njq2 ~= njq || niq2 ~= niq
    error('Grid and flux dimensions do not match!');
end

curlF = NaN * zeros(njq,niq);

[jc, ic] = deal(2:njq, 2:niq);

if weight_opt == 1
    dvdx = (Fy(jc,ic) - Fy(jc,ic-1)) ./ scqx(jc,ic);
    dudy = (Fx(jc,ic) - Fx(jc-1,ic)) ./ scqy(jc,ic);
elseif weight_opt == 2
    dvdx = 0; 
    dudy = 0;
end

curlF(jc,ic) = dvdx - dudy;