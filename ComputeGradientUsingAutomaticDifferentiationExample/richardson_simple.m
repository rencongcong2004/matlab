% 简单Richardson迭代，与chebyshev加速一起的
function [x, iter, res_history] = richardson_simple(A, b, x0, omega, tol, max_iter)
    x = x0;
    res_history = [];
    
    for iter = 1:max_iter
        r = b - A * x;
        res = norm(r);
        res_history = [res_history, res];
        
        if res < tol
            fprintf('Richardson迭代：%d次，残差=%.4e\n', iter, res);
            break;
        end
        
        x = x + omega * r;
    end
    
    if iter == max_iter
        fprintf('Richardson达到最大迭代：%d次，残差=%.4e\n', iter, res);
    end
end