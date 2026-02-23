n = 40;
for i=1:n
    x_values(i) = cos((2*i-1)*pi/(2*n));
end
fx = 1./(1+25.*x_values.^2);

y=ones(1,400);
y2=ones(1,400);
x1_values=linspace(-1, 1, 400);
for k=1:400
    y(k)=myBarycentricLagrange(x_values,fx,x1_values(k));
end
y2=1./(1+25.*x1_values.^2);
plot(x1_values,y2,'b');%蓝色曲线，是原函数图像
hold on
plot(x1_values,y,'r')  %红色曲线，是插值函数图像
error=y2-y;
Linf_error = norm(error, inf)



function [P_value] = myBarycentricLagrange(x_values,fx,x)
n=length(x_values);
    
% 计算wj
wj=ones(1,n);
for k=1:n
    indices=[1:k-1,k+1:n];
    wj(k)=1/prod(x_values(k)-x_values(indices));
end

%计算分母
Lk_denominators=ones(1,n);
for j=1:n
    Lk_denominators(j)=wj(j)/(x-x_values(j));
end
Wj=sum(Lk_denominators);

%计算分子
qj=ones(1,n);
for i=1:n
    qj(i)=Lk_denominators(i)*fx(i);
end
Qj=sum(qj);

P_value=Qj/Wj;
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
    y(k)=BarycentricLagrange(x_values,fx,x1_values(k));
end
y2=1./(1+25.*x1_values.^2);
plot(x1_values,y2,'b');%蓝色曲线，是原函数图像
hold on
plot(x1_values,y,'r')  %红色曲线，是插值函数图像
error=y2-y;
Linf_error = norm(error, inf)
%}