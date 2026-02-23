clear 

%误差
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
j1  = (1:n-1); %拟合节点个数加两端项等于多项式次数
xj1 = cos(pi.*j1./n)';

%等距节点
xj2 = linspace(-1,1,n-1)';

%前密后疏节点 
xj3 = linspace(0, 1, n-1)';  
l = 2;  
xj3 = 2 * xj3.^l - 1;  % 映射到 [-1, 1]

e1 = qde_VA(xj1,n,f,u,a,b,ahp,t);
e2 = qde_VA(xj2,n,f,u,a,b,ahp,t);
e3 = qde_VA(xj3,n,f,u,a,b,ahp,t);
fprintf('chebyshev节点:n = %d 时，最大误差为：%.4e\n', n, e1);
fprintf('等距节点:n = %d 时，最大误差为：%.4e\n', n, e2);
fprintf('前密后疏节点 :n = %d 时，最大误差为：%.4e\n', n, e3);