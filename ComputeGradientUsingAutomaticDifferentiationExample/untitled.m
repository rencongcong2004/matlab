function Q = LegendreBestpolyL2(x,m,f) %输入x为输出点的x值,m为输出多项式次数，h为内积（P，f）
%  生成Legendre多项式在输出点的值
x=x(:);
n = length(x);  
P=zeros(n,m); 
P(:,1) = 1; 
P(:,2) = x;
for k = 2:m-1
    P(:,k+1) = ((2*k-1)/k)*x.*P(:,k) - ((k-1)/k)*P(:,k-1); 
end

%计算C-C权重
n2=100;
c = zeros(n2+1,1);
c(1:2:end) = 2./(1-(0:2:n2).^2);
c(1) = c(1)/2;
c(end) = c(end)/2;
w = idct(c);
w = w';

% 生成Legendre多项式在C-C节点上的值
x2 =-cos(pi*(0:n2)'/n2); 
P2=zeros(n2+1,m); 
P2(:,1) = 1; 
P2(:,2) = x2;
for k = 2:m-1
    P2(:,k+1) = ((2*k-1)/k)*x2.*P2(:,k) - ((k-1)/k)*P2(:,k-1); 
end


h = zeros(m, 1);
G = zeros(m, 1);
for k = 1:m
    % 计算归一化系数
    h(k) = 2 / (2*(k-1) + 1);
    % 计算内积 (f, P_k)
    integrand = f(x2) .* P2(:, k);
    G(k) = (w * integrand) / h(k);    
end
    
    %  在输出点计算展开结果
Q = P * G;
end