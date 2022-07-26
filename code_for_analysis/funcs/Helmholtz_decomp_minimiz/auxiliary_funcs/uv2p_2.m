
function [u_p,v_p] = uv2p_2(u,v)
[u_p,v_p] = deal(zeros(size(u,1),size(v,2)));
u_p(:,2:end-1) = (u(:,1:end-1) + u(:,2:end)) / 2;
v_p(2:end-1,:) = (v(1:end-1,:) + v(2:end,:)) / 2;
end