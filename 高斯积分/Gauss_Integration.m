function  [G]=Gauss_Integration(n,rho,f,a,b)

P = cell(n+2, 1);
P{1} = @(x) 0;     % P_{-1} = 0（占位，不使用）
P{2} = @(x) 1;     % P_0 = 1

% 计算三项递推系数
for i = 2:n+1
    Pi = P{i};
    Pi_sq = @(x) Pi(x).^2;
    
    num_alpha = integral(@(x) rho(x) .* x .* Pi_sq(x), a, b);
    den_alpha = integral(@(x) rho(x) .* Pi_sq(x), a, b);
    
    if i == 2
        alpha_val = num_alpha / den_alpha;
        beta_val = integral(rho, a, b);  % β0
    else
        alpha_val = num_alpha / den_alpha;
        
        Pim1 = P{i-1};
        Pim1_sq = @(x) Pim1(x).^2;
        den_beta = integral(@(x) rho(x) .* Pim1_sq(x), a, b);
        beta_val = den_alpha / den_beta;
    end
    
    % 创建新的多项式函数（使用保存的系数）
    if i == 2
        % P_1(x) = (x - α0)P_0(x)
        P{i+1} = @(x) (x - alpha_val) .* Pi(x);
    else
        % P_k(x) = (x - α_{k-1})P_{k-1}(x) - β_{k-1}P_{k-2}(x)
        % 需要捕获当前的值
        current_alpha = alpha_val;
        current_beta = beta_val;
        current_Pi = Pi;          % P_{k-1}
        current_Pim1 = P{i-1};    % P_{k-2}

        P{i+1} = @(x) (x - current_alpha) .* current_Pi(x) - current_beta * current_Pim1(x);
        
    end
end


Pn_func = P{n+2};
syms x
% 使用Chebyshev节点采样（避免Runge现象）
num_points = 500;  % 增加采样点
cheb_nodes = cos(pi * (2*(1:num_points)-1) / (2*num_points));  % [-1,1]
x_vals = a + (b-a) * (cheb_nodes + 1) / 2;  % 映射到[a,b]

Pn_vals = zeros(1, num_points);
for i = 1:num_points
    Pn_vals(i) = Pn_func(x_vals(i));
end

% 使用稳健的拟合方法
poly_coeff = polyfit(x_vals, Pn_vals, n);
% Pn_sym = poly2sym(poly_coeff, x);

% 直接数值求根（避免符号运算误差）
roots_num = roots(poly_coeff);
roots_num = real(roots_num(abs(imag(roots_num)) < 1e-10));  % 只取实根
roots_num = roots_num(roots_num >= a & roots_num <= b);  % 在区间内的根
roots_num = sort(roots_num);




t=roots_num(:);
T = t.^(0:n-1); 
o=zeros(n, 1);
for u=1:n
    o(u)=integral(@(x) rho(x) .* x.^(u-1), a, b);
end
o=o(:);
d=(T'\o)';
fx=f(roots_num);
fx=fx(:);

G=d*fx
