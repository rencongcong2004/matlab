function Q = LegendreBestpolyL2(x,m,f) %x为输出点的x值,m为输出多项式次数
%  生成Legendre多项式在输出点的值
x=x(:);
n = length(x);  
P=zeros(n,m); 
P(:,1) = 1; 
P(:,2) = x;
for k = 2:m-1
    P(:,k+1) = ((2*k-1)./k).*x.*P(:,k) - ((k-1)./k).*P(:,k-1); 
end

%计算C-C权重
n2=100;
c = zeros(n2+1,1);%%%
c(1:2:n2+1) = 2./(1-(0:2:n2+1).^2);
w = imydct0(c);
w = w';

% 生成Legendre多项式在C-C节点上的值
x2 =-cos(pi*(0:n2)'/n2); 
P2=zeros(n2+1,m);
P2(:,1) = 1; 
P2(:,2) = x2;
for k = 2:m-1
     P2(:,k+1) = ((2*k-1)./k).*x2.*P2(:,k) - ((k-1)./k).*P2(:,k-1); 
end


h = zeros(m, 1);
G = zeros(m, 1);
for k = 1:m
    % 计算归一化系数
    h(k) = 2 / (2*(k-1) + 1);
    % 计算内积 (f, P_k)
    o = f(x2) .* P2(:, k);
    G(k) = (w * o) / h(k);    
end
    
    %  在输出点计算展开结果
Q = P * G;
end







%测试函数
% 
% clear; clc;
% 
% x = linspace(-0.9, 0.9, 100);
% 
% x = linspace(-1, 1, 10);
% 
% % 测试1：简单多项式（应该在低次数下精确重构）
% f1 = @(x) 1 + 2*x + 3*x.^2;  % 二次多项式
% m = 15;  
% Q1 = LegendreBestpolyL2(x, m, f1);
% 
% figure;
% subplot(2,2,1);
% plot(x, f1(x), 'b-', 'LineWidth', 2, 'DisplayName', 'f(x)=1+2x+3x^2');
% hold on;
% plot(x, Q1, 'ro', 'MarkerSize', 3, 'DisplayName', sprintf('Legendre展开 (m=%d)', m));
% legend; title('测试1：二次多项式');
% grid on;
% error1 = max(abs(f1(x) - Q1'));
% fprintf('测试1最大误差: %.10f\n', error1);
% 
% 
% %测试2：复杂函数
% f3 = @(x) exp(x);
% m3 = 15;
% Q3 = LegendreBestpolyL2(x, m3, f3);
% 
% subplot(2,2,3);
% plot(x, f3(x), 'b-', 'LineWidth', 2, 'DisplayName', 'f(x)=e^x');
% hold on;
% plot(x, Q3, 'ro', 'MarkerSize', 3, 'DisplayName', sprintf('Legendre展开 (m=%d)', m3));
% legend; title('测试3：指数函数');
% grid on;
% error3 = max(abs(f3(x) - Q3'));
% fprintf('测试3最大误差: %.10f\n', error3);
% 
% 
% % % 测试3：收敛性测试
% m_values = 4:2:16;
% errors = zeros(size(m_values));
% 
% subplot(2,2,4);
% for i = 1:length(m_values)
%     m = m_values(i);
%     Q = LegendreBestpolyL2(x, m, f3);
%     errors(i) = max(abs(f3(x) - Q'));
% end
% 
% semilogy(m_values, errors, 'o-', 'LineWidth', 2);
% xlabel('Legendre多项式次数 m');
% ylabel('最大误差');
% title('收敛性测试 (f(x)=e^x)');
% grid on;