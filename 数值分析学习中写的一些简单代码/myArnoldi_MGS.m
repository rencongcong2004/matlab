function[Q,H,normv] = myArnoldi_MGS(A,n,v)
m = size(A,1);
if nargin<3
    v=ones(m,1);
end


normv=norm(v);


if iscolumn(A)
    A=diag(A);
end


Q = zeros(m,n);
H = zeros(n+1,n);
v = v./norm(v); 
Q(:,1) = v;
for k = 1:n
    w = A*Q(:,k);
    for i = 1:k
        v = Q(:,i);
        H(i,k) = w'*v;
        w = w - H(i,k)*v;
    end
    H(k+1,k) = norm(w,2);
    Q(:,k+1) = w / H(k+1,k);
    if H(k+1,k) < 1e-12
        Q = Q(:,1:k); 
        H = H(1:k,1:k);
        break
    end
end



%测试函数

% m = 22; rng(0); A = rand (m);
% A = 5* speye (m) - (tril(A,3) - tril(A,-2));
% figure ; spy (A)
% r = ones(m,1); % -- test problem
% x0 = zeros(m,1); r0 = r - A*x0; e0 = norm(r0);
% n = 10; % choose n = 5, 10, 15, 20 and compare
% [Q,H,k] = myArnoldi_MGS(A,r0,n);
% bet = norm(r0); e1 = 0*H(:,1); e1(1) = bet;
% t = H \ e1;
% x1 = x0 + Q(:,1:k)*t;
% r1 = r - A*x1; e1 = norm(r1);
% fprintf('初始残差: %.6f\n', e0);
% fprintf('迭代次数: %d\n', k);
% fprintf('GMRES后残差: %.6f\n', e1);
% fprintf('改进倍数: %.2f倍\n', e0/e1);
