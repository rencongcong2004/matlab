clear


syms x
u = exp(sin(pi*x));
ahp = 0.1;
f =  -ahp.*diff(u, x, 2) + u;
f= matlabFunction(f);
u = matlabFunction(u);

a = 1;  b = 1;  ahp = 0.1;
t = linspace(-1,1,2e3)';%测试节点
n = 20;

%chebyshev节点
j1  = (0:n); %拟合节点个数加两端项等于多项式次数
xj = cos(pi.*j1./n)';


v0 = ones(n+1,1); 
xj = xj(:);
m = length(xj);
%这里需要调整G内积使得Q,d2Q的线性组合矩阵正交

Gprod = @(x,y) ( -ahp*x(2*m+1:3*m) + x(1:m) )' * ( -ahp*y(2*m+1:3*m) + y(1:m) );



[Q,~,d2Q,H] = zz_confArnoldi(xj, v0, n,Gprod);

A = -ahp*d2Q + Q;


orth_test = A' * A;
diag_vals = diag(orth_test);
A_norm = A * diag(1./sqrt(diag_vals));

cond(A_norm)
size(A_norm)

