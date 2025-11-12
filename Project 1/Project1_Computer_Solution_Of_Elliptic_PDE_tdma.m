clc;
clear all;

%% Solving the governing equation using line by line solver by making tdma structure

% Parametrs 
h = 250; % w/m^2-k
lambda = 5;    % w/m-k
Tinf = 573;    % k
L = 0.08;      % length in m 

grid_sizes = [11, 21, 41];
max_iter = 100000;
tolerance = 1e-6;

% S=making null arrays for results 
cpu_time = zeros(1, length(grid_sizes));
iterations = zeros(1, length(grid_sizes));
all_residuals = cell(1, length(grid_sizes));
all_T = cell(1, length(grid_sizes));

for g = 1:length(grid_sizes)
    N = grid_sizes(g);
    dx = L/(N-1);
    
    % Initialize temperature field initially
    T = zeros(N,N);
    
    % Apply boundary conditions
    T(1, 1:N) = 423;   %Top
    T(1:N, N) = 473;   % Right boundary  
    T(N, 1:N) = 323;   %Bottom
    
    % Left boundary convection
    for i = 2:N-1
        T(i,1) = (lambda/dx * T(i,2) + h * Tinf) / (lambda/dx + h);
    end
    
    tic;
    [T, residuals, iter] = line_by_line_solver(T, N, dx, h, lambda, Tinf, max_iter, tolerance);
    cpu_time(g) = toc;
    
    all_residuals{g} = residuals;
    iterations(g) = iter;
    all_T{g} = T;
    
    fprintf('Grid %dx%d: %d iterations, %.4f seconds\n', N, N, iter, cpu_time(g));
end

%% Plot temperature contour for finest grid

figure;
contourf(flipud(all_T{3}), 20);
colorbar;
title('Temperature Field - 41x41 Grid - Line by Line');
xlabel('Column Index');
ylabel('Row Index');

%% visualise the convergence heistory for different grid points

figure;
for grid_num = 1:length(grid_sizes)
    plot(1:length(all_residuals{grid_num}), all_residuals{grid_num}, 'LineWidth', 1.5);
    hold on;
end

legend('11x11', '21x21', '41x41 ');
xlabel('Number of Iterations ');
ylabel('Residual ');
title(' Convergence History For Different Grid Points');
grid on;

%% Plot CPU time vs grid points

figure;
plot([121, 441, 1681], cpu_time, 'o-', 'LineWidth', 2);
xlabel('Total Grid Points');
ylabel('CPU Time (s)');
grid on;
title('Line-by-Line Method - CPU Time vs Grid Points');
