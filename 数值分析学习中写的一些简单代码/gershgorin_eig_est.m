% Gershgorin 圆盘定理进行估计最大最小特征值
function [eig_lower_est, eig_upper_est] = gershgorin_eig_est(A)
% GERSHGORIN_EIG_EST 基于Gershgorin圆盘定理估计矩阵的最小/最大特征值
% 输入：
%   A - 待估计的n阶方阵（实对称/非对称均可）
% 输出：
%   eig_lower_est - 估计的最小特征值（实对称：特征值下界；一般矩阵：实部下界）
%   eig_upper_est - 估计的最大特征值（实对称：特征值上界；一般矩阵：实部上界）
% 说明：
%   1. 实对称矩阵：特征值为实数，估计值为严格的特征值上下界
%   2. 一般矩阵：特征值可能为复数，输出特征值实部的估计上下界

% 输入合法性检查
[n, m] = size(A);
if n ~= m
    error('输入必须为方阵！');
end

% 计算Gershgorin圆盘的中心（对角元）和半径（行非对角元绝对值和）
centers = diag(A);          % 圆盘中心（对角元）
radii = zeros(n, 1);        % 圆盘半径
for i = 1:n
    radii(i) = sum(abs(A(i, :))) - abs(centers(i)); % 每行非对角元绝对值和
end

% 估计特征值上下界
is_symmetric = isequal(A, A'); % 判断是否为实对称矩阵
if is_symmetric
    % 实对称矩阵：特征值为实数，直接取中心±半径的极值
    eig_lower_est = min(centers - radii);
    eig_upper_est = max(centers + radii);
else
    % 一般矩阵：输出特征值实部的上下界（复数特征值无直接"最大/最小"，聚焦实部）
    eig_lower_est = min(real(centers) - radii);
    eig_upper_est = max(real(centers) + radii);
end

% 仅输出估计的最小/最大特征值（无冗余打印，仅返回结果）
end





% 测试代码
% % 示例1：实对称矩阵（核心场景）
% A = [4, -1, 0; 
%      -1, 4, -1; 
%      0, -1, 4];
% [min_est, max_est] = gershgorin_eig_est(A);
% fprintf('实对称矩阵：\n');
% fprintf('估计最小特征值：%.4f\n', min_est);
% fprintf('估计最大特征值：%.4f\n', max_est);
% fprintf('真实最小特征值：%.4f\n', min(eig(A)));
% fprintf('真实最大特征值：%.4f\n', max(eig(A)));
% 
% % 示例2：一般非对称矩阵
% B = [2, -1, 3; 
%      1, 5, -2; 
%      -4, 2, 1];
% [min_est_B, max_est_B] = gershgorin_eig_est(B);
% fprintf('\n一般非对称矩阵（特征值实部）：\n');
% fprintf('估计最小实部：%.4f\n', min_est_B);
% fprintf('估计最大实部：%.4f\n', max_est_B);
% fprintf('真实特征值实部范围：[%.4f, %.4f]\n', min(real(eig(B))), max(real(eig(B))));