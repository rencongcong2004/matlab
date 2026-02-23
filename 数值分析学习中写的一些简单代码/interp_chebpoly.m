function p = interp_chebpoly(x,y,t) 
n = length(x);  
m = length(t); 

T = zeros(n,n); 
T(:,1) = 1; 
T(:,2) = x;
for k = 2:n-1
    T(:,k+1) = 2*x.*T(:,k) - T(:,k-1); 
end

cj = T\y; 

% cond(T), 
% j0 = (0:n-1).'; 
% T = cos(j0*j0'*pi/(n-1)); 
% cj = T\y; 




% cj = zz_idct2(y); 

T = zeros(m,n); 
T(:,1) = 1; 
T(:,2) = t;
for k = 2:n-1
    T(:,k+1) = 2*t.*T(:,k) - T(:,k-1); 
end

p = T*cj; 



end