%共轭梯度法
function [  epn, k ]=myCG(A,b,x0)
x0=x0(:);
b=b(:);
k = 0;
r0 = b - A*x0;
maxI = 10e3;
epn = 1;
Tol = 10e-10;
p0 = r0;
while k < maxI && epn > Tol
    d = A*p0;
    a = (r0' * p0)/(p0'*d); %或点积运算dot(r0, r0) ./ dot(r0, d);

    x1 = x0 + a*p0;
    r1 = r0 - a*d;

    f = -(r1'*d)/(p0'*d);
    p1 = r1 + f*p0;

    r0 = r1;
    p0 = p1;

    x0 = x1;
    k = k + 1;
    epn = norm(r0,inf);
end
fprintf (['\n CG method stop at itrstep=%d,' ...
    'and residual=%e\n'],k,epn);



% %测试代码
% clear
% clc
% A = [-4,1,1,1;1,-4,1,1;1,1,-4,1;1,1,1,-4];
% b = [1;1;1;1];
% x0 = [1;0;0;0];% x0 = [0;0;0;0];
% mygrad_descent(A,b,x0)
% myCG(A,b,x0)