clear
n = 100;
for i=1:n
    x_values(i) = cos((2*i-1)*pi/(2*n));
end
fx = 1./(1+25.*x_values.^2);

y=ones(1,400);
y2=ones(1,400);
x1_values=linspace(-1, 1, 400);
for k=1:400
    y(k)=myNewton_function(x_values,fx,x1_values(k));
end
y2=1./(1+25.*x1_values.^2);
plot(x1_values,y2,'b');%蓝色曲线，是原函数图像
hold on
plot(x1_values,y,'r')  %红色曲线，是插值函数图像
error=y2-y;
Linf_error = norm(error, inf)


function [P_value] = myNewton_function(x_values,fx,x)
n=length(x_values);
x_values=x_values(:);
fx=fx(:);
N=zeros(n,n);
N(:,1)=fx;
for j=2:n
    for i=j:n
        N(i,j)=(N(i,j-1)-N(i-1,j-1))/(x_values(i)-x_values(i-j+1));
    end
end

t=x-x_values;
p=ones(1,n);
p(1)=1;
for k=2:n
    p(k)=prod(t(1:k-1));
end

M=diag(N);
P_value=M'*p';
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
    y(k)=Newton_function(x_values,fx,x1_values(k));
end
y2=1./(1+25.*x1_values.^2);
plot(x1_values,y2,'b');%蓝色曲线，是原函数图像
hold on
plot(x1_values,y,'r')  %红色曲线，是插值函数图像
error=y2-y;
Linf_error = norm(error, inf)
%}