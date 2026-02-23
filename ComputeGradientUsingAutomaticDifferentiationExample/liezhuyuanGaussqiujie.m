function [x] =liezhuyuanGaussqiujie (A,b)
format long
n=size(A,1);
B=A;
P=eye(n);
L=eye(n);
U=eye(n);
M=eye(n);

for k=1:n-1
    [~,j]=max(abs(U(k:n,k)));
    j=k+j-1;
    if j~=k
        temp=A(k,:);
        A(k,:)=A(j,:);
        A(j,:)=temp;
        Mtemp=M(k,:);
        M(k,:)=M(j,:);
        M(j,:)=Mtemp;
    end
    A(k+1:n,k)=A(k+1:n,k)/A(k,k);
    A(k+1:n,k+1:n)=A(k+1:n,k+1:n)- A(k+1:n,k)* A(k,k+1:n);
    P=M*P;
end

for i=1:n-1
    L(i+1:n,i)=A(i+1:n,i);
end

for j=1:n
    U(j,j:n)=A(j,j:n);
end
h=P*b;
y=L\h;
x=U\y;
end


