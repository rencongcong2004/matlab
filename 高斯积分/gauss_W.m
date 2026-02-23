%gauss_chebpoly
function [xp,Ap] = gauss_W(n1,w,a,b)
%

if nargin<3
    a = -1; b = 1;
else
    xx = @(t) (b-a)/2*t + (a+b)/2;
    w = @(t) w( xx(t) ) ;
end



% step 1 compute orth poly (n degree) p = sum_k ck*Tk

tol = 1e-11; 

A = zeros(n1,n1); 
for k  = 1 : n1
    Tk = @(zz) get_chebpoly(zz,k-1);
    for j = 1 : k

        Tj =  @(zz) get_chebpoly(zz, j-1);

        Ij = integral(@(x) Tk(x).*Tj(x).*w(x), -1, 1 );
        % Ij = quad(@(x) Tk(x).*Tj(x).*w(x), -1, 1 , tol);

        A(k,j) = Ij ;
    end
end

A = A + tril(A,-1).';

% kA = cond(A),

rhs =  zeros(n1,1);

Tn1 = @(zz) get_chebpoly(zz, n1);


for k = 1:n1
    Tk = @(zz) get_chebpoly(zz, k-1);

    rhs(k) = integral(@(x) Tk(x).*Tn1(x).*w(x), -1,1 );
end



ck = A \  -rhs(:) ;
% step 2 compute roots of p
xp = roots_chebpoly([ ck; 1]) ;


% step 3 compute quad coeffs
rhs = zeros(n1,1);
for k = 1:n1
    Tk = @(zz) get_chebpoly(zz, k-1);

    A(:,k) = Tk(xp) ;
    rhs(k) = integral(@(zz) w(zz).*Tk(zz), -1,1);
end

Ap = (A.') \ rhs(:); % gauss coef


xp = (b-a)*xp/2 + (a+b)/2;
Ap = (b-a)/2.*Ap;


end

% % % % ==================subroutine 1 =============================

function T = get_chebpoly(x, n)

T1 = 1;  T2= x;

for j = 2:n
    T3 = 2*x.*T2- T1;
    T1 = T2; T2 = T3;
end

if n == 0
    T = T1;
elseif n==1
    T = T2;
else
    T = T3;
end
end
% % % % ==================subroutine 2 =============================


function [zz,x] = roots_chebpoly(ck)
%     colleague matrix of chebyshev series
%     find roots of p = sum ck*Tk(x) in [-1,1] // k = 1->n
%     referrence :  Yuji, et al, 2016, math. comp.
% ON THE STABILITY OF COMPUTING POLYNOMIAL ROOTS VIA CONFEDERATE LINEARIZATIONS
% II = find(abs(ck) > 8e-16,1,'last');
% if isempty(II)
%     error('Input coeffs dk are all zeros');
% elseif II>1
%     ck = ck(1:II);
% end
ck = ck/norm(ck);
n = length(ck)-1;
% cheb_hess(n) ;
% Hessenberg matrix for chebyshev polys
% size = (n+1,n)
alp = 1/2*ones(n+1,1);
H =   diag(alp(1:n),1) + diag(alp(2:n+1),-1)  ;
H(2,1) = 1;
H = H(:,1:n);
meth = 2;

if n > 500
    H = sparse(H);
    meth = 3;
    warning(['roots for %d-th chebyshev series may costly and not' ...
        'efficient '], n);
end


switch meth
    case 1
        %         may ill-conditon if d(n+1) is small
        ck = ck./ck(n+1) ;
        H1 = H(1:n,1:n);
        cc = ck(1:n).*H(n+1,n);
        H1(:,n) = H1(:,n) - cc;
        zz = eig(H1);

    case 2
        % general eig : //   S^-1*A = H1
        A = H(1:n,1:n);
        S = eye(n);         S(n,n) = ck(n+1);
        A(:,n) =  A(:,n)*ck(n+1) - H(n+1,n)*ck(1:n);
        [~,x] = eig(A, S, 'qz' );
        zz = sort(diag(x));
    case 3
        % meth 2; by eigs sparse
        A = H;  A(end,:) = [];
        S = speye(n);         S(n,n) = ck(n+1);
        A(:,n) =  A(:,n)*ck(n+1) - H(n+1,n)*ck(1:n);
        x = eigs(A, S, n );
        zz = sort(x);
end
% select real roots
% zz,
% zz  = zz(abs(zz)<=1);
% zz = zz(abs( imag(zz) ) < 2e-15);
% zz = sort(real(zz));

end



 