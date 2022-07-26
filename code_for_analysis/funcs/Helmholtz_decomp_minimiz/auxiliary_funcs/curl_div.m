function [curl,div] = curl_div(u,v,sc)

[scqx,scqy] = deal(sc.scqx, sc.scqy); 
[scpx,scpy] = deal(sc.scpx, sc.scpy); 
[scux,scuy] = deal(sc.scux, sc.scuy); 
[scvx,scvy] = deal(sc.scvx, sc.scvy); 

% [cqx,cqy] = deal(cxy.cqx, cxy.cqy); 
% [cpx,cpy] = deal(cxy.cpx, cxy.cpy); 
% [cux,cuy] = deal(cxy.cux, cxy.cuy); 
% [cvx,cvy] = deal(cxy.cvx, cxy.cvy);

curl = NaN * zeros(size(u));
% JDM-1, IDM-1
dvdx = (v(:,2:end) - v(:,1:end-1)) ./ scqx(2:end-1,2:end-1);
dudy = (u(2:end,:) - u(1:end-1,:)) ./ scqy(2:end-1,2:end-1);
% curl (q-)
curl(2:end,2:end) = dvdx - dudy;

%  div (p-)
div = calc_div_HYCOM(u,v,scvx,scuy,scpx.*scpy,1);

end