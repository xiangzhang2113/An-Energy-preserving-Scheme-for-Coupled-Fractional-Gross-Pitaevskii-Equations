%Table 1. Numerical errors and convergence orders in time with h = 0.01 at final time T = 0.5. 
h = 1/100;
tau = [1/40, 1/80, 1/160, 1/320, 1/640];
alpha_u = [2, 2, 2, 1.8, 1.2];
alpha_v = [2, 1.8, 1.5, 2, 1.8];
P = length(alpha_u);
K = length(tau);

% 用 cell 数组保存所有结果
rate_l2_all = cell(P,1);
rate_max_all = cell(P,1);

for p = 1:P
    fprintf('Processing alpha pair #%d: alpha_u = %.2f, alpha_v = %.2f\n', p, alpha_u(p), alpha_v(p));
    
    % 初始化 cell 保存该组 alpha 的结果
    u_II_cell = cell(K,1);
    epsilon_cell = cell(K-1,1);
    maxepsilon_cell = cell(K-1,1);

    % 计算数值解
    [u_II_cell{1}, v_II_cell{1}] = fgp_solver(alpha_u(p), alpha_v(p), tau(1), h);
    [u_II_cell{K}, v_II_cell{K}] = fgp_solver(alpha_u(p), alpha_v(p), tau(K), h);
    for n = 2:K-1
        [u_II_cell{n}, v_II_cell{n}] = fgp_solver(alpha_u(p), alpha_v(p), tau(n), h);
        epsilon_cell{n-1} = sqrt(h * sum( abs(u_II_cell{n-1} - u_II_cell{n}).^2 ));
        maxepsilon_cell{n-1} = max(abs(u_II_cell{n-1} - u_II_cell{n}));
    end
    epsilon_cell{K-1} = sqrt(h * sum( abs(u_II_cell{K-1} - u_II_cell{K}).^2 ));
    maxepsilon_cell{K-1} = max(abs(u_II_cell{K-1} - u_II_cell{K}));

    % 计算收敛阶
    rate_l2 = zeros(K-2,1);
    rate_max = zeros(K-2,1);
    for n = 1:K-2
        rate_l2(n) = log2(epsilon_cell{n} / epsilon_cell{n+1});
        rate_max(n) = log2(maxepsilon_cell{n} / maxepsilon_cell{n+1});
    end

    % 存储结果
    rate_l2_all{p} = rate_l2;
    rate_max_all{p} = rate_max;

    % 可选：显示本轮结果
    disp(['L2 convergence rate for p = ', num2str(p), ':']);
    disp(rate_l2);
    disp(['Max error convergence rate for p = ', num2str(p), ':']);
    disp(rate_max);
end








