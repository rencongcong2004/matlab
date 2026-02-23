%qde问题各种节点误差

clear 


syms x
u = exp(sin(pi*x));
ahp = 0.1;
f =  -ahp.*diff(u, x, 2) + u;
f= matlabFunction(f);
u = matlabFunction(u);


a = 1;  b = 1;  ahp = 0.1;
t = linspace(-1,1,2e3)';%测试节点
uj = u(t);

n = 20;
%chebyshev节点 效果最好
j1  = (1:n-1); %拟合节点个数加两端项等于多项式次数
xj1 = cos(pi.*j1./n);


%等距节点
xj2 = linspace(-1,1,n-1)';


%前密后疏节点 
xj3 = linspace(0, 1, n-1)';  
l = 2;  
xj3 = 2 * xj3.^l - 1;  % 映射到 [-1, 1]



ck1 = qde_chebyshev(n,xj1,ahp,f,a,b);
ck2 = qde_chebyshev(n,xj2,ahp,f,a,b);
ck3 = qde_chebyshev(n,xj3,ahp,f,a,b);

[T,~,~] = mydiff_chebpolyT(t,n) ;

Un1 = T*ck1;
Un2 = T*ck2;
Un3 = T*ck3;



e1 = norm(Un1-uj,inf);
e2 = norm(Un2-uj,inf);
e3 = norm(Un3-uj,inf);

fprintf('chebyshev节点:n = %d 时，最大误差为：%.4e\n', n, e1);
fprintf('等距节点:n = %d 时，最大误差为：%.4e\n', n, e2);
fprintf('前密后疏节点 :n = %d 时，最大误差为：%.4e\n', n, e3);



% %画图
% figure;
% 
% subplot(1,2,1);
% plot(xj, Un, 'r--', 'LineWidth', 2, 'DisplayName', '谱方法解');  hold on;
% plot(xj,uj,'b-', 'LineWidth', 2, 'DisplayName', '精确解');
% grid on;
% 
% subplot(1,2,2);
% error = Un - uj;
% plot(xj, error, 'r-', 'LineWidth', 2);
% grid on;
% 




% % %画图验证前密后疏节点 在x轴上位置
% figure;
% plot(x_points, zeros(size(x_points)), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% xlim([ -1.1*max(x_points), 1.1*max(x_points)]);
% ylim([-0.1, 0.1]);
% xlabel('x 值');
% title('节点在区间上的位置'); 