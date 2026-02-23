clear 
tic
%误差曲线
syms x
u = exp(sin(pi*x));
ahp = 0.1;
f =  -ahp.*diff(u, x, 2) + u;
f= matlabFunction(f);
u = matlabFunction(u);

a = 1;  b = 1;  ahp = 0.1;
t = linspace(-1,1,2e3)';%测试节点
m = 40;
E1 = [];
E2 = [];
E3 = [];

for i = 2:m
    n = i;

    %chebyshev节点 效果最好
    j1  = (1:n-1); %拟合节点个数加两端项等于多项式次数
    xj1 = cos(pi.*j1./n);
    
    
    %等距节点
    xj2 = linspace(-1,1,n-1)';
    
    
    %前密后疏节点 
    xj3 = linspace(0, 1, n-1)';  
    l = 2;  
    xj3 = 2 * xj3.^l - 1;  % 映射到 [-1, 1]

    e1 = qde_VA(xj1,n,f,u,a,b,ahp,t);
    e2 = qde_VA(xj2,n,f,u,a,b,ahp,t);
    e3 = qde_VA(xj3,n,f,u,a,b,ahp,t);
    E1 = [E1,e1];
    E2 = [E2,e2];
    E3 = [E3,e3];
end

n = linspace(2,m,m-1)';


figure;
plot(n, E1, 'r-', 'DisplayName', 'Chebyshev节点误差', 'LineWidth', 2);
hold on;
plot(n, E2, 'b-', 'DisplayName', '等距节点误差', 'LineWidth', 2);
plot(n, E3, 'g-', 'DisplayName', '前密后疏节点误差', 'LineWidth', 2);

% 添加标签和标题
ylim([0, 0.1]); % 只显示 0-0.1 的范围
xlabel('多项式次数 n', 'FontSize', 12);
ylabel('最大误差', 'FontSize', 12);
title('不同节点类型的误差比较', 'FontSize', 14);
grid on;

% 显示图例
legend('Location', 'best');

toc