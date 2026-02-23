function [xj,wj] = zz_Riemann_Stieltjes(g, n, xj)
%
% compute : S f(x) d g(x)  by chebpts interp
% if isempty(g)
%     g = @(x)  x;
% end

% n = 5;
% g = @(x)  x;
if nargin<3 | isempty(xj)
    xj =  cos( (0:n)'*pi/n);
    % xj = cos( ((0:n)'+0.5)*pi/(n+1) );
end

r = zeros(n+1,1);

% dT = diff( chebpoly(0:n) ) ;

for k = 1:n+1
    %
    % mj = @(zz) zz.^(k-1).*g(zz);
    % mj = T(:,k);
    % mj2 = integral(@(zz) mj(zz).*dg(zz), a, b);
    % r(k ) =  mj2;

    bd = g(1) - (-1).^(k-1)*g(-1);

    dTk = @(zz) diff_chebpoly(zz,k-1);
    % dTk = dT(:,k);

    inn = integral(@(zz) g(zz).*dTk(zz), -1, 1 );


    r(k ) =  bd - inn ;

end

T = zz_chebpoly(xj,n);
wj = T' \ r ;

% wj = imydct0(r);


[xj,ia] = sort(xj);
wj = wj(ia);


end


% % % % ==================subroutine 1 =============================

function dT = diff_chebpoly(x, n)

T1 = 1;         T2 = x;
dT1 = x*0 ;  dT2=  1+x*0;

for j = 2:n
    T3 = 2*x.*T2 - T1;
    dT3 = 2*x.*dT2 + 2*T2 - dT1;
    T1 = T2;      T2 = T3;
    dT1 = dT2;  dT2 = dT3;
end

if n == 0
    dT = dT1;
elseif n==1
    dT = dT2;
else
    dT = dT3;
end
end

% % % % ==================subroutine 2 =============================
function y = imydct0(x)  % ==>    inverse of mydct0
% C = cos(jk*pi/(m-1)), j,k = 0 -> m-1
% y = C \ x
[m,~] = size(x);  x([1,m],:) = 0.5*x([1,m], :);
% y = mydct0(x);
x(m+1:2*m-2,:) = 0;
y = fft(x);  y = real(y(1:m,:));
%
y([1,m],:) = 0.5*y([1,m], :);  y = y*(2/(m-1));
end
