clear 
%qde问题各种节点误差曲线
tic

syms x
u = exp(sin(pi*x));
ahp = 0.1;
f =  -ahp.*diff(u, x, 2) + u;
f= matlabFunction(f);
u = matlabFunction(u);


a = 1;  b = 1;  ahp = 0.1;
t = linspace(-1,1,2e3)';%测试节点
uj = u(t);
E1 = [];
E2 = [];
E3 = [];

k = 2e3;
m = 40;
for i = 2:m
    n = i;
    %chebyshev节点 效果最好
    j1  = (1:n-1); %拟合节点个数加两端项等于多项式次数
    xj1 = cos(pi.*j1./n);
    
    
    %等距节点
    xj2 = linspace(-1,1,k)';
    
    
    %前密后疏节点 
    xj3 = linspace(0, 1,k)';  
    l = 2;  
    xj3 = 2 * xj3.^l - 1;  % 映射到 [-1, 1]
    
    
    
    ck1 = qde_chebyshev(n,xj1,ahp,f,a,b);
    ck2 = qde_chebyshev(n,xj2,ahp,f,a,b);
    ck3 = qde_chebyshev(n,xj3,ahp,f,a,b);
    
    [T,~,~] = mydiff_chebpolyT(t,n) ;
    
    Un1 = T*ck1;
    Un2 = T*ck2;
    Un3 = T*ck3;
    
    
    
    
    e1 = norm(Un1-uj,inf);
    e2 = norm(Un2-uj,inf);
    e3 = norm(Un3-uj,inf);

    E1 = [E1,e1];
    E2 = [E2,e2];
    E3 = [E3,e3];
end

n = linspace(2,m,m-1)';


figure;
semilogy(n, E1, 'r-', 'DisplayName', 'Chebyshev节点误差', 'LineWidth', 2);
hold on;
semilogy(n, E2, 'b-', 'DisplayName', '等距节点误差', 'LineWidth', 2);
semilogy(n, E3, 'g-', 'DisplayName', '前密后疏节点误差', 'LineWidth', 2);

% 添加标签和标题
ylim([0, 0.1]); % 只显示 0-0.1 的范围
xlabel('多项式次数 n', 'FontSize', 12);
ylabel('最大误差', 'FontSize', 12);
title('不同节点类型的误差比较', 'FontSize', 14);
grid on;

% 显示图例
legend('Location', 'best');

toc