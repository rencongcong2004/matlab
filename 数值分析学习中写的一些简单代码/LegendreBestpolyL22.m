function Q = LegendreBestpolyL2(x,m,f) % x为节点值,m为输出多项式次数，f为要逼近函数

%  生成Legendre多项式在节点处的值
x=x(:);
n = length(x);  
% P = zeros(n,m); 
% P(:,1) = 1; 
% P(:,2) = x;
% for k = 2:m-1
%     % P(:,k+1) = ((2*k+1)/(k+1))*x.*P(:,k) - (k/(k+1))*P(:,k-1); 
% end

P = zz_legpoly(x,m-1); 

% %计算C-C权重
% n2 = 100;
% c = zeros(n2+1,1);
% c(2:2:n2+1) = 2./(1-(2:2:n2+1).^2);
% w = imydct2(c);
% w = w';
% 
% % 生成Legendre多项式在C-C节点上的值
% x2 = cos(pi*(0:n2)'/n2); 

n2 = 100; 
[x2,w] = zz_chebpts(n2); 

P2 = zz_legpoly(x2,m-1); 

% P2=zeros(n2+1,m);
% P2(:,1) = 1; 
% P2(:,2) = x2;
% for k = 2:m-1
%     P2(:,k+1) = ((2*k-1)/(k ))*x2.*P2(:,k) - ( (k-1)/(k ))*P2(:,k-1);
% end


h = zeros(m, 1);
G = zeros(m, 1);
for k = 1:m
    % 计算归一化系数
    h(k) = 2 / (2*(k-1) + 1);
    % 计算内积 (f, P_k)
    o = f(x2) .* P2(:, k);
    G(k) = (w * o) / h(k);    
end
    
    %  在输出点计算展开结果
Q = P * G;
end



% % % % % % % % % % % % subroutine=======================


function [V,dV,d2V] = zz_legpoly(x,m,a,b)
% for given points x, and intger m
% it comput P(x) , P'(x), P"(x) on x, with degree from 0\to m
x = x(:); N = length(x);
if (nargin < 2),   m = N-1;   end


if (nargin < 3),   a = -1; b = 1;    end


w = 2/(b-a);
z = (x - (b+a)/2).*w;
V = zeros(N,m+1);   dV = V;   d2V = V;
V(:,1) = 1;  V(:,2) = z;    dV(:,2) = w;

switch nargout 
    case 1
    for j = 2:m
        alp =  (2*j-1)./j;   bet = (j-1)./j;
        V(:,j+1) = alp.*z.*V(:,j) - bet.*V(:,j-1) ;
    end
    
    case 2
    for j = 2:m
        alp =  (2*j-1)./j;   bet = (j-1)./j;
        V(:,j+1) = alp.*z.*V(:,j) - bet.*V(:,j-1) ;
        dV(:,j+1) =  alp.*(w.*V(:,j) +  z.*dV(:,j) ) - bet.*dV(:,j-1);
    end

    case 3
    for j = 2:m
        alp =  (2*j-1)./j;   bet = (j-1)./j;
        V(:,j+1) = alp.*z.*V(:,j) - bet.*V(:,j-1) ;
        dV(:,j+1) =  alp.*(w.*V(:,j) +  z.*dV(:,j) ) - bet.*dV(:,j-1);
        d2V(:,j+1) =  alp.*(2*w.*dV(:,j) +  z.*d2V(:,j) ) - bet.*d2V(:,j-1);
    end
end

end



% % % % % % % % % % % % subroutine=======================


% %计算C-C权重
% n2 = 100;
% c = zeros(n2+1,1);
% c(2:2:n2+1) = 2./(1-(2:2:n2+1).^2);
% w = imydct2(c);
% w = w';
% 
% % 生成Legendre多项式在C-C节点上的值
% x2 = cos(pi*(0:n2)'/n2); 

function [x,w,If] = zz_chebpts(n,a,b,f)
%   Clenshaw-Curtis for If = int_{a}^{b} f(x) dx
%   Input:
%               n # of quad pts
%               f // funciton handle
%   Output:
%               x: chebyshev-lobatto pts on [a,b]
%               w: weight of integral
%
x = -cos(pi*(0:n).'/n);         w = 0*x;  j0 = (0:2:n).';       % (n+1) chebpts
w(1:2:n+1) = 2./(1-j0.^2);  w = imydct0(w);  w = w.';

%  % trans [-1,1] to [a,b]
if (nargin > 1)
    x = (b-a)*x/2 + (a+b)/2;
    w = (b-a)/2.*w;
end

if nargin < 4
    If = [];
    return
else
    If = w*f(x);
end

end



% -------------------------------------------------------
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

