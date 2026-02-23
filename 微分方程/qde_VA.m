function [e] = qde_VA(xj,n,f,u,a,b,ahp,t)
v0 = ones(n-1,1); 
xj = xj(:);
%这里需要调整G内积使得Q,d2Q的线性组合矩阵正交
[Q,~,d2Q,H] = zz_confArnoldi(xj, v0, n);

A = -ahp*d2Q + Q;

v_norm = norm(v0);  
[u1] = zz_confArnoldi_inv(H, -1, v_norm) ;
[u2] = zz_confArnoldi_inv(H, 1, v_norm) ;
U = [u1;u2];

A = [A;U];

%计算右端项
fj = f(xj);
fj = fj(:);
fj = [fj;a;b];
c = A\fj;

uj = u(t);
[Q ] = zz_confArnoldi_inv(H, t, v_norm) ;

Un1 = Q*c;
e = norm(Un1-uj,inf);