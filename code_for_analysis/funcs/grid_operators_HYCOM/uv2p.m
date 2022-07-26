function [u_p,v_p] = uv2p(u,v)

[u_p,v_p] = deal(NaN * zeros(size(u)));
u_p(:,1:end-1) = (u(:,1:end-1) + u(:,2:end)) / 2;
v_p(1:end-1,:) = (v(1:end-1,:) + v(2:end,:)) / 2;

end