%简化的Chebyshev实现（基于Richardson迭代）
function [x, iter, res_history] = simple_chebyshev_iter(A, b, x0, lambda_min, lambda_max, tol, max_iter)
    % 简化的Chebyshev实现（基于Richardson迭代）
    
    % 计算最优参数
    omega = 2 / (lambda_max + lambda_min);
    rho = (lambda_max - lambda_min) / (lambda_max + lambda_min);
    d = rho^2;
    
    % 初始化
    x = x0;
    r = b - A * x;
    res = norm(r);
    res_history = [res];
    iter = 0;
    
    if res < tol
        return;
    end
    
    % 第一步：普通Richardson
    x_old = x;
    x = x + omega * r;
    r = b - A * x;
    res = norm(r);
    res_history = [res_history, res];
    iter = 1;
    
    % Chebyshev加速
    for k = 2:max_iter
        % 计算Chebyshev参数
        if k == 2
            beta_k = 1 / (2 * (1 - d));
        else
            beta_k = 1 / (1 - d/4);
        end
        
        % 更新解
        x_new = beta_k * (x + omega * r) + (1 - beta_k) * x_old;
        
        % 计算残差
        r = b - A * x_new;
        res = norm(r);
        res_history = [res_history, res];
        
        % 更新变量
        x_old = x;
        x = x_new;
        iter = k;
        
        % 收敛检查
        if res < tol
            break;
        end
    end
    
    fprintf('Chebyshev迭代：%d次，残差=%.4e\n', iter, res);
end






% 测试，cheby加速与CG和最速下降法的比较
% clear
% clc
% 
% % 
% % n = 100; 
% % beta=2.1;
% % 
% % 
% % main_diag = beta*ones(n,1);% 主对角向量
% % sub_diag = 1*ones(n-1,1); % 上/下对角向量
% % % 构造三对角矩阵
% % A = diag(main_diag) + diag(sub_diag,1) + diag(sub_diag,-1);
% % 
% % b=rand(n,1);
% % x0 = rand(n,1);
% % 
% % mygrad_descent(A,b,x0)
% % myCG(A,b,x0)
% % 
% %
% clear; clc; close all;
% 
% %% 1. 构造三对角矩阵
% n = 5000; 
% beta = 3;
% main_diag = beta * ones(n,1);
% sub_diag = 1 * ones(n-1,1);
% A = diag(main_diag) + diag(sub_diag,1) + diag(sub_diag,-1);
% 
% b = rand(n,1);
% x0 = rand(n,1);
% x_exact = A \ b;
% 
% %% 2. 计算A的特征值范围
% lambda_min_A = beta + 2 * cos(n * pi / (n + 1));  % λ_min ≈ 0.1002
% lambda_max_A = beta + 2 * cos(pi / (n + 1));      % λ_max ≈ 4.0998
% 
% fprintf('矩阵A特征值范围：[%.6f, %.6f]\n', lambda_min_A, lambda_max_A);
% fprintf('条件数 κ = %.4f\n', lambda_max_A/lambda_min_A);
% 
% %% 3. 验证Chebyshev实现
% fprintf('\n=== 验证Chebyshev加速 ===\n');
% 
% % 最优松弛因子
% omega = 2 / (lambda_max_A + lambda_min_A);
% fprintf('最优松弛因子 ω = %.6f\n', omega);
% 
% % 运行Chebyshev
% tol = 1e-9;
% max_iter = 10000;
% 
% % 方法1：使用MATLAB的chebyshev半迭代（如果有的话）
% % 方法2：手动实现简化版
% 
% [x_cheb, iter_cheb, res_cheb] = simple_chebyshev_iter(A, b, x0, lambda_min_A, lambda_max_A, tol, max_iter);
% 
% %% 4. 与其他方法对比
% fprintf('\n=== 普通Richardson迭代 ===\n');
% [x_rich, iter_rich, res_rich] = richardson_simple(A, b, x0, omega, tol, max_iter);
% 
% fprintf('\n=== 共轭梯度法(CG) ===\n');
% [x_cg, flag_cg, relres_cg, iter_cg, resvec_cg] = pcg(A, b, tol, max_iter, [], [], x0);
% 
% %% 5. 绘图
% figure('Position', [100, 100, 900, 400]);
% 
% subplot(1,2,1);
% semilogy(res_rich, 'b-', 'LineWidth', 2); hold on;
% semilogy(res_cheb, 'r--', 'LineWidth', 2);
% semilogy(resvec_cg, 'g:', 'LineWidth', 2);
% xlabel('迭代次数');
% ylabel('残差');
% legend('Richardson', 'Chebyshev', 'CG', 'Location', 'best');
% title('收敛历史');
% grid on;
% 
% subplot(1,2,2);
% bar([iter_rich, iter_cheb, iter_cg]);
% set(gca, 'XTickLabel', {'Richardson', 'Chebyshev', 'CG'});
% ylabel('迭代次数');
% title('收敛所需迭代次数');
% grid on;
% 
% %% 辅助函数

