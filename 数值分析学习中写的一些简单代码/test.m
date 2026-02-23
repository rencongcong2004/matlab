clear,
LW = 'LineWidth'; 
% test fun
f = @(x) 1./(1 + 25.*x.^2); 
a = 0;  b = 1; 

m = 2e3; 
t = linspace(a,b,m).'; 

rho = @(x) sqrt(1 - x.^2); 

W = sparse(1:m,1:m,rho(t)); 

% % test lsq approximation 
n1 = 32 ; 
ft = f(t); 

% case 1 Vandermonde type basis 
T = t.^(0:n1); 
% cj = T\ft; 
% A = T'*W*T; 
% RHS = T'*W*ft; 
% cj = A\RHS; 
[QQ,RR] = qr(sqrt(W)*T, 0); 

cj = RR\(QQ'*sqrt(W)*ft); 
 
cond(RR), 

pt = T*cj; 

e1 = norm(pt-ft,inf), 

% % case 2  chebyshev
T = mychebpolyT(t,n1); 
% T = chebpoly(0:n1 );
% T = T(t);
[QQ,RR] = qr( sqrt(W)*T, 0 ); 

cj = RR\(QQ'*( sqrt(W)*ft ) ) ; 

pt2 = T*cj; 

cond(RR),
% 
e2 = norm(pt2 - ft,inf), 
% 
% case 3  V + A
v = sqrt(rho(t)); 
[Q,H,normv] = myArnoldi(t, n1,v); 
%cond(Q)
d = Q'*( v.*ft ) ; 

d = d./normv; 

pt3 = mypolyvalA(d,H,t); 

e3 = norm(pt3 - ft,inf), 
 
% % 



