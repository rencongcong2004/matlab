n=19;
% a = -1; b = 0.3;
a=6;b=8;

rho= @(x) 1./sqrt((1-x.^2));
% rho =  @(x) 1./sqrt((1-x.^4));  %这两个rho无明显差异，
% rho = @(x) x.^8;        %这个差异较大


f = @(x) 1./(1 + 25.*x.^2); 
G1=integral(@(x) rho(x) .* f(x), a, b);

tic
[xp1,Ap1]= gauss_Legendre(n,rho,a,b);
ft1=f(xp1);
ft1=ft1(:);
Ap1=Ap1(:)';
G2=Ap1*ft1;
e1=norm(G1-G2,inf), 
toc


tic
[xp2,Ap2]=  gauss_W(n,rho,a,b);
ft2=f(xp2);
ft2=ft2(:);
Ap2=Ap2(:)';
G3=Ap2*ft2;

e2=norm(G1-G3,inf), 
toc