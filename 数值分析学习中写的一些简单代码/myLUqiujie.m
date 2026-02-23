function [x] = myLUqiujie(A,b)
format long
U=A;
n=size(A,1);
L=eye(n);
for k=1:n-1
    L(k+1:n,k)=U(k+1:n,k)/U(k,k);
    U(k+1:n,k+1:n)=U(k+1:n,k+1:n)- L(k+1:n,k)* U(k,k+1:n);
    U(k+1:n,k)=0;
end

y=L\b;
x=U\y;

