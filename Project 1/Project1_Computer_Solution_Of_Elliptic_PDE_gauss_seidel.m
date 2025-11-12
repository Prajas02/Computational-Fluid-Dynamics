clc;
clear all;
%% Solving the governing equation using Gauss-Seidel iterative approach

% Params
h = 250; % W/m^2K
lambda = 5; % W/m-K
Tinf = 573; % K
L = 0.08; % length of domain in m
grid_sizes = [11, 21, 41];
max_iter = 100000;
tolerance = 1e-6;

% for comparison between different grid sizes make a loopp and solve for
% each different grid points step by step

% S=making null arrays for results 
cpu_time = zeros(1, length(grid_sizes));
iterations = zeros(1, length(grid_sizes));
all_residuals = cell(1, length(grid_sizes));
all_T = cell(1, length(grid_sizes));

for grid_num = 1:length(grid_sizes)
    N = grid_sizes(grid_num);
    dx = L/(N-1);
    
    % Initialize temperature field
    T =  zeros(N,N);
    
    % Apply boundary conditions
    T(1, 1:N) = 423;      % Top boundary
    T(1:N, N) = 473;      % Right boundary
    T(N, 1:N) = 323;      % Bottom boundary
    
    % Initial guess for left boundary (convective)
    for i = 2:N-1
        T(i,1) = (lambda * T(i,2) + h * dx * Tinf) / (lambda + h * dx);
    end
    
    % start noting down the cpu time
    tic;
    [T, residuals, iter] = gauss_seidel_method(T, N, dx, h, lambda, Tinf, max_iter, tolerance);
    cpu_time(grid_num) = toc;
    
    all_residuals{grid_num} = residuals;
    iterations(grid_num) = iter;
    all_T{grid_num} = T;
    
    fprintf('Grid %dx%d: %d iterations, %.4f seconds\n', N, N, iter, cpu_time(grid_num));
end

%% Plot temperature contour for finest grid

figure;
contourf(flipud(all_T{3}), 40);
colorbar;
title('Temperature Field - 41x41 Grid');
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
plot([121, 441, 1681], cpu_time, 'o-', 'LineWidth', 1.5 );
xlabel(' Total number of Grid Points ');
ylabel('CPU Time in seconds' );
title('CPU Time vs Total Number of Grid Points ');
grid on;