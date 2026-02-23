function  T = mychebpolyT(x,m) 
x = x(:);
n = length(x);  
T = zeros(n,m+1); 
T(:,1) = 1; 
T(:,2) = x;
for k = 2:m
    T(:,k+1) = 2*x.*T(:,k) - T(:,k-1); 
end