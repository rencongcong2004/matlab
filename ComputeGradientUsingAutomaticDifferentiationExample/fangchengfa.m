n = 40;
for i=1:n
    x_values(i) = cos((2*i-1)*pi/(2*n));
end
fx = 1./(1+25.*x_values.^2);

y=ones(1,n);
y2=ones(1,n);
x1_values=ones(1,n);
for k=1:n
    x=1/k;
    x1_values(k)=1/k;
    y(k)=myfangchengfa(x_values,fx,x);
    y2=1./(1+25.*x1_values.^2);
end
plot(x1_values,y2,'b');
hold on
plot(x1_values,y,'r')
error=y2-y;
Linf_error = norm(error, inf)



function [P_value] = myfangchengfa(x_values,fx,x)
n=length(x_values);
x_values=x_values(:);
b=fx(:);
A=zeros(n);
A(1:n,1)=1;
for k=2:n
    A(1:n,k)=x_values.^(k-1);
end

P=ones(1,n);
a=liezhuyuanGaussqiujie (A,b);
for i=1:n
    P(i)=x^(i-1)*a(i);
end
P_value=sum(P);

end


%{
n = 40;
for i=1:n
    x_values(i) = cos((2*i-1)*pi/(2*n));
end
fx = 1./(1+25.*x_values.^2);

y=ones(1,n);
y2=ones(1,n);
x1_values=ones(1,n);
for k=1:n
    x=1/k;
    x1_values(k)=1/k;
    y(k)=fangchengfa(x_values,fx,x);
    y2=1./(1+25.*x1_values.^2);
end
plot(x1_values,y2,'b');
hold on
plot(x1_values,y,'r')
error=y2-y;
Linf_error = norm(error, inf)
%}




%{
n = 40;
for i=1:n
    x_values(i) = cos((2*i-1)*pi/(2*n));
end
fx = 1./(1+25.*x_values.^2);

y=ones(1,n);
x1_values=linspace(-1, 1, n);
y2=1./(1+25.*x1_values.^2);
for k=1:n
    y(k)=fangchengfa(x_values,fx,x1_values(k));
end
plot(x1_values,y2,'b');
hold on
plot(x1_values,y,'r')
error=y2-y;
Linf_error = norm(error, inf)
%}