function  f=mypolyvalA(d,H,x,norm_g0)
d = d(:);
x = x(:);
n = length(d);
m = length(x);
P=zeros(m,n);
v = ones(m,1)./norm_g0;
P(:,1)=v;

for j=1:n-1
    w=zeros(m,1);
    for i=1:j
        w=w+H(i,j)*P(:,i);
    end
    T=x.*P(:,j)-w;
    P(:,j+1)=T./H(j+1,j);
end

f=P*d;
