clear


syms x
u = exp(sin(pi*x));
ahp = 0.1;
f =  -ahp.*diff(u, x, 2) + u;
f= matlabFunction(f);
u = matlabFunction(u);

a = 1;  b = 1;  ahp = 0.1;
t = linspace(-1,1,2e3)';%测试节点
E1 = [];

 k=60;
for i = 2:k
    n = i;
    
    %chebyshev节点
    j1  = (1:n-1); %拟合节点个数加两端项等于多项式次数
    xj = cos(pi.*j1./n)';
    
    
    v0 = ones(n-1,1); 
    xj = xj(:);
    m = length(xj);
    %这里需要调整G内积使得Q,d2Q的线性组合矩阵正交
    
    Gprod = @(x, y) ( -ahp*x(2*m+1:3*m) + x(1:m) )' * ...
                    ( -ahp*y(2*m+1:3*m) + y(1:m) ) + ...
                    x(3*m+1:end)' * y(3*m+1:end);  % 边界条件部分
    
    
    
    [Q,~,d2Q,QQ,H,norm_v0] = zz_confArnoldi2(xj, v0, n,Gprod);
    size(Q);
    A = -ahp*d2Q + Q;
    
    A = [A;QQ];
    cond(A);

  
    % v_norm = sqrt(norm(v0));  
    fj = f(xj);
    fj = fj(:);
    fj = [fj;a;b];
    c = A\fj;
    
    uj = u(t);
    [Q ] = zz_confArnoldi_inv(H, t,  norm_v0) ;
    
    Un1 = Q*c;
    e = norm(Un1-uj,inf);
    E1 = [E1,e];
end
% 
% 
n = linspace(2,k,k-1)';
figure;

semilogy(n,E1, 'r-', 'DisplayName', 'Chebyshev节点误差', 'LineWidth', 2);
hold on;
% 添加标签和标题
ylim([0, 0.1]); % 只显示 0-0.1 的范围
xlabel('多项式次数 n', 'FontSize', 12);
ylabel('最大误差', 'FontSize', 12);
title('不同节点类型的误差比较', 'FontSize', 14);
grid on;

% 显示图例
legend('Location', 'best');

