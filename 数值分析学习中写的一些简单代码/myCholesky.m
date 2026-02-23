function [] = myCholesky(A)
format long
B=A;
n=size(A,1);
L=eye(n);
for k=1:n
    A(k,k)=sqrt(A(k,k));
    A(k+1:n,k)= A(k+1:n,k)/A(k,k);
    for j=k+1:n
        A(j:n,j)=A(j:n,j)-A(j:n,k)*A(j,k);
    end
end

for i=1:n
    L(i:n,i)=A(i:n,i);
end

U=L*L';
D=U-B;
disp(L)
disp(U)
disp(D)


%A=[4,-2,4,2;-2,10,-2,-7;4,-2,8,4;2,-7,4,7];b=[8;2;16;6];
%y=L\b;
%x=L'\y;
%disp(x)