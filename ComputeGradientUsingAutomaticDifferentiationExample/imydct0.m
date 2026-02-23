function y = imydct0(x)
m = size(x,1); 
x([1,m],:) = 0.5*x([1,m], :);
y = mydct0(x); 
y([1,m],:) = 0.5*y([1,m], :); 
y = y*(2/(m-1));
end

function y = mydct0(x)
[m,n] = size(x); 
y = zeros(2*(m-1),n); 
y(1:m,:) = x(1:m,:);
y = fft(y); 
y = real(y(1:m,:));
end