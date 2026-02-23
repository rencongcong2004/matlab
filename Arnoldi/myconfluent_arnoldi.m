function [G,H,norm_g0] = myconfluent_arnoldi(xj,m,v)
%xj为拟合节点，m是多项式次数
n = length(xj);
xj = xj(:);

G = zeros(n,m);
dG = zeros(n,m);
H = zeros(m+1,m);
%初始化
g0 = v;
norm_g0 = sqrt(g0' * g0);
G(:, 1) = g0 / norm_g0;
dG(:, 1) = zeros(n, 1);

x = diag(xj);

for i = 1:m
    w1 = x*G(:,i);
    w2 = G(:,i) + x*dG(:,i);
    for k = 1:i
        f = [G(:,k); dG(:,k)];
        g = [w1; w2];
        H(k,i) = f'*g;


        v1 = H(k,i).*G(:,k);
        v2 = H(k,i).*dG(:,k);

        w1 = w1 - v1;
        w2 = w2 - v2;

    end

     H(i+1,i) =  sqrt(w1' * w1 + w2' * w2);

    
    G(:,i+1) = w1 ./ H(i+1,i);
    dG(:,i+1) = w2 ./ H(i+1,i);
    if H(i+1,i) < 1e-12
        G = G(:,1:k);
        H = H(1:k,1:k);
        break
    end
    
end

G = [G;dG];