function [c]=qde_chebyshev(n,xj,ahp,f,a,b)

%n是多项式次数，xj是拟合节点，ahp，a，b是参数，f是右端项函数


%除两端项之外的系数矩阵
[T,~,dT2] = mydiff_chebpolyT(xj,n) ;
H = -ahp.*dT2 + T;


%计算两端项
x1 = 1;
x1 = x1(:);

[V1,~,~] = mydiff_chebpolyT(-1.*x1,n) ;
[V2,~,~] = mydiff_chebpolyT(x1,n) ;


%计算最终系数矩阵

A = [H;V1;V2];

%计算右端项
fj = f(xj);
fj = fj(:);
fj = [fj;a;b];

%求解
c = A\fj;