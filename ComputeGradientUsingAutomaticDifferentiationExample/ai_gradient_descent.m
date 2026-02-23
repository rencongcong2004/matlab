%梯度下降法
function [x, fval, k, history] = ai_gradient_descent(f, grad_f, x0, alpha, max_iter, tol, verbose)
% 通用梯度下降法优化器
% 输入：
%   f: 目标函数句柄，f(x) 返回标量
%   grad_f: 梯度函数句柄，grad_f(x) 返回梯度向量
%   x0: 初始点
%   alpha: 步长策略 - 标量（固定步长）或函数句柄 alpha(k)
%   max_iter: 最大迭代次数（默认：1000）
%   tol: 收敛容差（默认：1e-6）
%   verbose: 是否显示进度（默认：1，每100次显示）
% 输出：
%   x: 最优解
%   fval: 最优函数值
%   k: 实际迭代次数
%   history: 迭代历史记录 [k, fval, norm_grad, step_size]

    % ========== 参数处理与验证 ==========
    if nargin < 4
        error('至少需要4个参数：f, grad_f, x0, alpha');
    end
    
    if nargin < 5 || isempty(max_iter), max_iter = 1000; end
    if nargin < 6 || isempty(tol), tol = 1e-6; end
    if nargin < 7 || isempty(verbose), verbose = 1; end
    
    % 输入验证
    x0 = x0(:);  % 确保列向量
    if ~isa(f, 'function_handle')
        error('f必须是函数句柄');
    end
    if ~isa(grad_f, 'function_handle')
        error('grad_f必须是函数句柄');
    end
    
    % ========== 初始化 ==========
    x = x0;
    fval = f(x);                % 当前函数值
    grad = grad_f(x);           % 当前梯度
    norm_grad = norm(grad);     % 梯度范数
    k = 0;                      % 迭代计数器
    converged = false;          % 收敛标志
    conv_reason = '达到最大迭代次数'; % 默认收敛原因
    
    % 存储迭代历史（预分配内存提高效率）
    history = zeros(min(max_iter, 10000), 4);
    history_idx = 1;
    
    % 记录初始状态
    history(history_idx, :) = [k, fval, norm_grad, 0];
    history_idx = history_idx + 1;
    
    % ========== 显示初始信息 ==========
    if verbose >= 1
        fprintf('\n=== 梯度下降法开始 ===\n');
        fprintf('维度: %d\n', length(x0));
        fprintf('初始函数值: %.6e\n', fval);
        fprintf('初始梯度范数: %.6e\n', norm_grad);
        fprintf('最大迭代次数: %d\n', max_iter);
        fprintf('收敛容差: %.1e\n\n', tol);
        
        if verbose >= 2
            fprintf('迭代  函数值       梯度范数      步长\n');
            fprintf('----------------------------------------\n');
        end
    end
    
    % ========== 主迭代循环 ==========
    for k = 1:max_iter
        % 1. 获取步长
        if isa(alpha, 'function_handle')
            alpha_k = alpha(k);          % 自适应步长
        else
            alpha_k = alpha;             % 固定步长
        end
        
        % 步长安全检查
        if alpha_k <= 0 || ~isfinite(alpha_k)
            warning('步长无效: alpha_k = %.2e，使用默认值0.01', alpha_k);
            alpha_k = 0.01;
        end
        
        % 2. 更新解
        x_new = x - alpha_k * grad;
        
        % 3. 计算新状态
        fval_new = f(x_new);
        grad_new = grad_f(x_new);
        norm_grad_new = norm(grad_new);
        
        % 4. 检查是否发散
        if ~isfinite(fval_new) || ~isfinite(norm_grad_new)
            if verbose >= 1
                fprintf('警告：函数值或梯度非有限，迭代终止\n');
            end
            conv_reason = '函数值或梯度非有限';
            converged = false;
            break;
        end
        
        % 5. 检查Armijo条件（确保充分下降）- 可选
        % 如果函数值不下降，减小步长
        if fval_new > fval
            if verbose >= 2
                fprintf('迭代 %d: 函数值增加，考虑减小步长\n', k);
            end
            % 这里可以实现回溯线搜索
        end
        
        % 6. 记录历史
        history(history_idx, :) = [k, fval_new, norm_grad_new, alpha_k];
        history_idx = history_idx + 1;
        
        % 7. 显示进度
        if verbose >= 2 && (mod(k, 100) == 0 || k <= 10)
            fprintf('%4d  %.3e  %.3e  %.2e\n', ...
                k, fval_new, norm_grad_new, alpha_k);
        end
        
        % 8. 收敛性检查（多重条件）
        x_change = norm(x_new - x);
        f_change = abs(fval_new - fval);
        grad_norm = norm_grad_new;
        
        % 收敛条件：梯度足够小 或 参数变化小 或 函数值变化小
        if grad_norm < tol                     % 梯度收敛
            conv_reason = '梯度范数达到容差';
            converged = true;
            break;
        elseif x_change < tol * (1 + norm(x))  % 参数收敛
            conv_reason = '参数变化达到容差';
            converged = true;
            break;
        elseif f_change < tol * (1 + abs(fval)) % 函数值收敛
            conv_reason = '函数值变化达到容差';
            converged = true;
            break;
        end
        
        % 9. 为下一次迭代准备
        x = x_new;
        fval = fval_new;
        grad = grad_new;
        norm_grad = norm_grad_new;
    end
    
    % ========== 后处理与输出 ==========
    % 截断历史记录
    history = history(1:history_idx-1, :);
    
    % 最终输出
    if verbose >= 1
        fprintf('\n=== 迭代结束 ===\n');
        fprintf('最终迭代次数: %d\n', k);
        fprintf('最终函数值: %.6e\n', fval);
        fprintf('最终梯度范数: %.6e\n', norm_grad);
        
        if converged
            fprintf('✓ 收敛原因: %s\n', conv_reason);
        else
            fprintf('⚠ %s\n', conv_reason);
        end
        fprintf('总计算时间: 需要实际添加计时\n');
    end
    
    % 绘制收敛曲线（可选）
    if verbose >= 1 && ~isempty(history)
        figure;
        subplot(2,2,1);
        semilogy(history(:,1), history(:,2), 'b-', 'LineWidth', 1.5);
        xlabel('迭代次数'); ylabel('函数值(对数)'); title('函数值收敛');
        grid on;
        
        subplot(2,2,2);
        semilogy(history(:,1), history(:,3), 'r-', 'LineWidth', 1.5);
        xlabel('迭代次数'); ylabel('梯度范数(对数)'); title('梯度范数收敛');
        grid on;
        
        subplot(2,2,3);
        plot(history(:,1), history(:,4), 'g-', 'LineWidth', 1.5);
        xlabel('迭代次数'); ylabel('步长'); title('步长变化');
        grid on;
        
        subplot(2,2,4);
        plot(history(2:end,1), diff(history(:,2)), 'm-', 'LineWidth', 1.5);
        xlabel('迭代次数'); ylabel('函数值变化'); title('函数值下降量');
        grid on;
    end
end


%测试代码
% clear
% clc
% % % 1定义Rosenbrock函数（香蕉函数）
% rosenbrock = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
% grad_rosenbrock = @(x) [ -400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
%                           200*(x(2)-x(1)^2) ];
% 
% x0 = [-1.2; 1];
% 
% % 固定步长
% alpha = 0.001;  
% 
% [x_opt, fval, k, history] = ai_gradient_descent(...
%     rosenbrock, grad_rosenbrock, x0, alpha, 5000, 1e-6, 2);
% % % % 衰减步长策略
% % alpha_func = @(k) 0.1 / (1 + 0.01*k);  % 衰减步长
% % 
% % [x_opt, fval, k, history] = ai_gradient_descent(...
% %     rosenbrock, grad_rosenbrock, x0, alpha_func, 5000, 1e-6, 1);




% % 2生成线性回归数据
% m = 100; n = 5;
% X = randn(m, n);
% w_true = randn(n, 1);
% y = X * w_true + 0.1 * randn(m, 1);
% 
% % 定义最小二乘损失和梯度
% loss = @(w) 0.5 * norm(X*w - y)^2;
% grad_loss = @(w) X' * (X*w - y);
% 
% w0 = zeros(n, 1);
% alpha = 0.01;  % 固定步长
% 
% [w_opt, fval, k] = ai_gradient_descent(loss, grad_loss, w0, alpha, 1000, 1e-8, 1);
% 
% % 与解析解比较
% w_exact = X \ y;
% fprintf('与解析解误差: %e\n', norm(w_opt - w_exact));