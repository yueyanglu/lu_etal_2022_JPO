function [u,v] = p2uv(u_p,v_p)
% 
% Re-define u and v velocity fields on p-grid onto u-/v- grids.
% 
[u,v] = deal(NaN * zeros(size(u_p)));
u(:,2:end) = (u_p(:,1:end-1) + u_p(:,2:end)) / 2;
v(2:end,:) = (v_p(1:end-1,:) + v_p(2:end,:)) / 2;

end