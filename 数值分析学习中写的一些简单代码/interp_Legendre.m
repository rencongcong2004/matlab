function p = interp_Legendre(x,y,t) 
n = length(x);  
m = length(t); 

P=zeros(n,n); 
P(:,1) = 1; 
P(:,2) = x;
for k = 2:n-1
    P(:,k+1) = ((2*k+1)/(k+1))*x.*P(:,k) - (k/(k+1))*P(:,k-1); 
end

cj = P\y; 

P=zeros(m,n); 
P(:,1) = 1; 
P(:,2) = t;
for k = 2:n-1
    P(:,k+1) = ((2*k+1)/(k+1))*t.*P(:,k) - (k/(k+1))*P(:,k-1); 
end

p=P*cj;
end


%{
测试代码
clear
n = 822;
i0 = (1:n).'; 
x_values = cos((2*i0-1)*pi/(2*n));
fx = 1./(1+25.*x_values.^2);

tic,
m= 4e3;
x1_values=linspace(-1, 1, m)';
y =  interp_Legendre(x_values,fx,x1_values);
y2 = 1./(1+25.*x1_values.^2);

toc,

plot(x1_values,y2,'b');%蓝色曲线，是原函数图像
hold on
plot(x1_values,y,'r')  %红色曲线，是插值函数图像
error=y2-y;
Linf_error = norm(error, inf)
%}