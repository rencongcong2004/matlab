clear
n = 822;
i0 = (1:n).'; 
x_values = cos((2*i0-1)*pi/(2*n));
fx = 1./(1+25.*x_values.^2);

tic,
m= 4e3;
x1_values=linspace(-1, 1, m)';
y =mycreate_lagrange_function(x_values,fx,x1_values);
y2=1./(1+25.*x1_values.^2);
toc,

plot(x1_values,y2,'b');%蓝色曲线，是原函数图像
hold on
plot(x1_values,y,'r')  %红色曲线，是插值函数图像
error=y2-y;
Linf_error = norm(error, inf)



function [P_value] =mycreate_lagrange_function(x_values,fx,x)
n=length(x_values);
Lk_denominators=ones(n,1);
m = length(x);
fenzi = ones(m,n); ell = fenzi; 
for k = 1:n
    indices=[1:k-1, k+1:n];
    % 计算分母
    Lk_denominators(k)=prod(x_values(k)-x_values(indices));
    % 计算分子
    fenzi(:,k) = prod(x-x_values(indices).', 2);
    ell(:,k) = fenzi(:,k)./Lk_denominators(k); 
end
    % P_value = sum(fx.*fenzi./Lk_denominators);
    P_value = ell*fx;

end


%{
n = 40;
for i=1:n
    x_values(i) = cos((2*i-1)*pi/(2*n));
end
fx = 1./(1+25.*x_values.^2);

y=ones(1,400);
y2=ones(1,400);
x1_values=linspace(-1, 1, 400);
for k=1:400
    y(k)=create_lagrange_function(x_values,fx,x1_values(k));
end
y2=1./(1+25.*x1_values.^2);

plot(x1_values,y2,'b');%蓝色曲线，是原函数图像
hold on
plot(x1_values,y,'r')  %红色曲线，是插值函数图像
error=y2-y;
Linf_error = norm(error, inf)
%}
