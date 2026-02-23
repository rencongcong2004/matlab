function  [T,dT1,dT2] = mydiff_chebpolyT(x,m) 
x=x(:);
n = length(x);  
T = zeros(n,m+1); 
T(:,1) = 1; 
T(:,2) = x;
for k = 2:m
    T(:,k+1) = 2.*x.*T(:,k) - T(:,k-1); 
end

dT1 = zeros(n,m+1);
dT1(:,1) = 0; 
dT1(:,2) = 1;
for k = 2:m
    dT1(:,k+1) = 2*x.*dT1(:,k) - dT1(:,k-1) + 2.*T(:,k); 
end


dT2 = zeros(n,m+1);
dT2(:,1) = 0; 
dT2(:,2) = 0;
dT2(:,3) = 4;
for k = 3:m
    dT2(:,k+1) = 2*x.*dT2(:,k) - dT2(:,k-1) + 4.*dT1(:,k); 
end