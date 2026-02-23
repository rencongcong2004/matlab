function [Q,dQ,d2Q,QQ,H,norm_v0 ] = zz_confArnoldi2(X, v, n, varargin)
%
% this code compute bases of Krylov subspace [e,X*e,X^2*e,...,X^n*e]
% span(Q) = span(1,x,...,x^n);
% span(dQ) =  d/dx( span(1,x,...,x^n)  ) ;
% span(d2Q) = d2/dx2 ( span(1,x,...,x^n);
%
%
% references:
%
%  [1] Q. Niu, H. Zhang and Y. Zhou, Confluent Vandermonde with Arnoldi,
%         Appl. Math. Lett., 135 (2023), 108420.
%  [2] Lei-Hong Zhang, Ya-Nan Zhang, Linyi Yang, Yifu Wu,  Multivariate confluent
%        Vandermonde with G-Arnoldi and applications, 2024

mj = abs(X(end) +  X(1)) ;
if mj > 0.1
    Xmean = mean(X);
    X = X-Xmean;
else
    Xmean = 0;
end



m = length(X);


if isempty(v)
    v = ones(m,1);
end

% -----------------

if nargin<4
    Gproduct = @(x,y) x'*y;
else
    Gproduct = varargin{:};
end
% example 1: hermite interpolant

% example 2:  2p ode    u + a*u"  = f


Q = zeros(m,n);
dQ = Q;  d2Q = Q;
H = zeros(n+1,n);
Q(:,1) = v;
% 构造完整向量 [Q; dQ; d2Q; 边界条件]
f = [Q(:, 1); dQ(:, 1); d2Q(:, 1); 1; 1];
mj = Gproduct(f,f);

norm_v0 = sqrt(mj);

f = f./norm_v0;
Q(:,1) = f(1:m,1);
dQ(:,1) = f(1+m:2*m,1);
d2Q(:,1) = f(2*m+1:3*m,1);
QQ = zeros(2,n);
QQ(1,1) = f(3*m+1,1);
QQ(2,1) = f(3*m+2,1);



for k = 1:n
    w = X.*Q(:,k);
    dw = Q(:,k) + X.*dQ(:,k);
    d2w = 2*dQ(:,k) + X.*d2Q(:,k);

    Qw1 = -1*QQ(1,k);
    Qw2 = QQ(2,k);

    %            modified Gram-Schmidt
    for i = 1:k
        %         v = Q(:,i);
        %         H(i,k) = v'*w;
        f = [Q(:,i); dQ(:,i); d2Q(:,i);QQ(1,i);QQ(2,i)];
        g = [w; dw; d2w;Qw1;Qw2];


        H(i,k) = Gproduct(f,g);


        w = w - Q(:,i)*H(i,k);
        dw = dw - dQ(:,i)*H(i,k);
        d2w = d2w - d2Q(:,i)*H(i,k);
        Qw1 = Qw1 - QQ(1,i)*H(i,k);
        Qw2 = Qw2 - QQ(2,i)*H(i,k);
    end

    %   re-orth

    % g = [w; dw; d2w];
    % for i = 1:k
    %     mj = Gproduct(f,g);
    %     w = w - Q(:,i) * mj;
    %     dw = dw - dQ(:,i) * mj;
    %     d2w = d2w - d2Q(:,i) * mj;
    %     H(i,k) = H(i,k) + mj;
    % end



    %     H(k+1,k) = norm(w,2);
    f = [w; dw; d2w;Qw1;Qw2];
    mj = Gproduct(f,f );
    H(k+1,k) = sqrt(mj);

    Q(:,k+1) = w / H(k+1,k);

    dQ(:,k+1) = dw / H(k+1,k);
    d2Q(:,k+1) = d2w / H(k+1,k);

    QQ(1,k+1) = Qw1 / H(k+1,k);
    QQ(2,k+1) = Qw2 / H(k+1,k);


    if H(k+1,k) < 2e-15
        fprintf('\n arnoldi process break down at step %d \n',k);
        H = H(1:k,1:k);
        Q = Q(:,1:k);
        dQ = dQ(:,1:k);
        d2Q = d2Q(:,1:k);
        break
    end

end

if Xmean~=0
    H = H + Xmean*eye(size(H));
end

end





