%gauss_Legendre
function [xp,Ap] = gauss_Legendre(n1,w,a,b)
%

if nargin<3
    a = -1; b = 1;
else
    xx = @(t) (b-a)/2*t + (a+b)/2;
    w = @(t) w( xx(t) ) ;
end

%计算系数矩阵A
A = zeros(n1,n1); 
for k  = 1 : n1
    Tk = @(zz) get_Legendre(zz,k-1);
    for j = 1 : k

        Tj =  @(zz) get_Legendre(zz, j-1);

        Ij = integral(@(x) Tk(x).*Tj(x).*w(x), -1, 1 );
        % Ij = quad(@(x) Tk(x).*Tj(x).*w(x), -1, 1 , tol);

        A(k,j) = Ij ;
    end
end

A = A + tril(A,-1).';

% kA = cond(A),

rhs =  zeros(n1,1);

Tn1 = @(zz) get_Legendre(zz, n1);


for k = 1:n1
    Tk = @(zz) get_Legendre(zz, k-1);

    rhs(k) = integral(@(x) Tk(x).*Tn1(x).*w(x), -1,1 );
end



ck = A \  -rhs(:) ;
% step 2 compute roots of p
xp = roots_Legendre([ ck; 1]) ;


% step 3 compute quad coeffs
rhs = zeros(n1,1);
for k = 1:n1
    Tk = @(zz) get_Legendre(zz, k-1);

    A(:,k) = Tk(xp) ;
    rhs(k) = integral(@(zz) w(zz).*Tk(zz), -1,1);
end

Ap = (A.') \ rhs(:); % gauss coef


xp = (b-a)*xp/2 + (a+b)/2;
Ap = (b-a)/2.*Ap;


end





function [zz] = roots_Legendre(ck)

ck = ck/norm(ck);
n = length(ck)-1;

j0 = (1:n); 
a = (2.*j0-1)./j0;
c=(j0-1)./j0;

a1=1./a;
a1=a1(1:n-1);

gama = c./a;
gama = gama(2:n);

H = diag(a1,1) + diag (gama,-1);

ck = ck./ck(n+1) ;
ck = ck';
u = (2*n-1)/n;
ck = ck./u;
H(n,:) = H(n,:)-ck(1:n);
zz = eig(H);
end




function T = get_Legendre(x, n)

T1 = 1;  T2= x;

for j = 2:n
    T3 = (2*j-1)*x.*T2- (j-1).*T1;
    T3 = T3/j;
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