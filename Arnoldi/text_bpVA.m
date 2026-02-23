clear
syms x
f =  1./(1+25.*x.^2);
df1 = diff(f, x);
f= matlabFunction(f);
df1= matlabFunction(df1);
m = 2000;
a = -1;  b = 1; 

n = 44; %多项式次数
m = 2e3;%拟合节点个数

%主程序
%等距节点
x_points = linspace(a,b,m).';%chebpoly节点上进行拟合

%右端项
fj = [f(x_points);df1(x_points)];

%系数矩阵
% Gprod = @(x,y) x(1:2*m)'*y(1:2*m); 
% v0 = ones(m,1); 
% [G,dG,~,H] = zz_confArnoldi(x_points, v0, n,Gprod);
% G = [G;dG];
v0 = ones(m,1); 
[G,H] =myconfluent_arnoldi(x_points, n,v0);


cond(G)
c = G'*fj;

%测试误差
t = linspace(-1,1,2*m)';%测试节点

ft = f(t);
v_norm = norm(v0); 

g = mypolyvalA(c,H,t,v_norm);

% [W,dW ] = zz_confArnoldi_inv(H,  t, v_norm) ; 
% g = W*c;


e = norm(g - ft,inf),


