
syms x
f =  1./(1+25.*x.^2);
df1 = diff(f, x);
f= matlabFunction(f);
df1= matlabFunction(df1);



m = 44; %多项式次数
n = 2e3;%拟合节点个数

% %主程序
% %chebshev节点
% j  = (0:n); 
% x_points = cos(pi.*j./n)';%chebpoly节点上进行拟合

%等距节点
x_points = linspace(-1,1,n)';


fj = f(x_points);
df1j = df1(x_points);

[T,dT1,~] = mydiff_chebpolyT(x_points,m) ;

A = [T;dT1];%系数矩阵
Fj = [fj;df1j];%右端项

c1 = A\Fj;

%计算误差
xj = linspace(-1,1,2*n)';%测试节点
[T2,~,~] = mydiff_chebpolyT(xj,m) ;

e2 = norm(T2*c1 - f(xj),inf);
fprintf('m = %d 次多项式，拟合点个数为m = %d 个时，最大误差为：%.4e\n', m, n,e2);






% %等距节点
% x_points = linspace(-1,1,n)';

% %前密后疏节点 该形式对这个问题效果较好
% xj1 = linspace(0, 1, n)';  
% l = 2;  
% x_points = 2 * xj1.^l - 1;  % 映射到 [-1, 1]
% %画图验证节点在x轴上位置
% figure;
% plot(x_points, zeros(size(x_points)), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% xlim([ -1.1*max(x_points), 1.1*max(x_points)]);
% ylim([-0.1, 0.1]);
% xlabel('x 值');
% title('节点在区间上的位置'); 