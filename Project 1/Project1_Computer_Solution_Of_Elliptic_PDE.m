clc;
clear all;

%% Solving the governing equation using Gaussian Elimination

% Initializing the boundary conditions
h = 250;       % W/m^2K
lambda = 5;    % W/m-K
Tinf = 573;    % K
L = 0.08;      % length of domain in m

grid_sizes = [11, 21, 41];

% Boundary condition
T_bottom = 323;
T_top = 423;
T_right = 473;

% making function for calculating teh time for different grid resolutions
% to solve everything at once

% S=making zero arrays for results 
cpu_time = zeros(1, length(grid_sizes));
iterations = zeros(1, length(grid_sizes));
all_residuals = cell(1, length(grid_sizes));
all_T = cell(1, length(grid_sizes));

for grid_num = 1:length(grid_sizes)
    N = grid_sizes(grid_num);
    
    % Start cpu time
    tic;
    
    % coefficient matrix and rhs_mat
    [A, b] = compute_coefficient_rhs_mat(N, T_bottom, T_top, T_right, h, lambda, Tinf);
    
    % Solve using Gaussian elimination
    interior_solution = gaussian_elimination(A, b);
    
    % Reconstruct full temperature field
    T = compute_temperature_field(interior_solution, N, T_bottom, T_top, T_right, h, lambda, Tinf);
    
    cpu_time(grid_num) = toc;
    all_T{grid_num} = T;
    
    fprintf('Grid %dx%d: %.4f seconds\n', N, N, cpu_time(grid_num));
end

%% Plot temperature contour for finest grid
figure;
contourf(all_T{3}, 40);
colorbar;
title('Temperature Field - 41x41 Grid (Gaussian Elimination)');
xlabel('X');
ylabel('Y');
axis equal tight;

%% Plot CPU time vs grid points
figure;
plot([121, 441, 1681], cpu_time, 'o-', 'LineWidth', 2);
xlabel('Total Number of Grid Points');
ylabel('CPU Time (seconds)');
title('Gaussian Elimination: CPU Time vs Grid Size');
grid on;

%% Function to compute coefficient matrix and RHS vector

function [A, b] = compute_coefficient_rhs_mat(N, T_bottom, T_top, T_right, h, lambda, Tinf)
    % Number of interior points
    interior_points = (N-2) * (N-2);
    dx = 0.08/(N-1);
    
    A = zeros(interior_points, interior_points);
    b = zeros(interior_points, 1);
    
    for j = 2:N-1     % y direction
        for i = 2:N-1 % x 
            
            index = (j-2)*(N-2) + (i-1); % done in class by (i,j) = 1
            
            % Center
            A(index, index) = -4;
            
            % Left neighbour
            if i > 2
                left_index = (j-2)*(N-2) + (i-2);
                A(index, left_index) = 1;
            else
                % Left boundary convective
                A(index, index) = A(index, index) + 1/(1 + (h * dx / lambda));
                b(index) = b(index) - ((h * dx / lambda) * Tinf) / (1 + (h * dx / lambda));
            end
            
            % Right neighbour
            if i < N-1
                right_index = (j-2)*(N-2) + i;
                A(index, right_index) = 1;
            else
                % Right boundary
                b(index) = b(index) - T_right;
            end
            
            % Bottom neighbour
            if j > 2
                bottom_index = (j-3)*(N-2) + (i-1);
                A(index, bottom_index) = 1;
            else
                % Bottom boundary
                b(index) = b(index) - T_bottom;
            end
            
            % Top neighbour
            if j < N-1
                top_index = (j-1)*(N-2) + (i-1);
                A(index, top_index) = 1;
            else
                % Top boundary
                b(index) = b(index) - T_top;
            end
        end
    end
end

%% Function to compute temperature field

function T = compute_temperature_field(interior_solution, N, T_bottom, T_top, T_right, h, lambda, Tinf)
    T = zeros(N, N);
    dx = 0.08/(N-1);
    
    % solve for interior points from solution vector
    for j = 2:N-1
        for i = 2:N-1
            index = (j-2)*(N-2) + (i-1);
            T(j, i) = interior_solution(index); 
        end
    end
    
    % Apply boundary conditions
    T(1, :) = T_bottom;
    T(N, :) = T_top;    
    T(:, N) = T_right;   
    
    % Compute left boundary using Neumann condition
    for j = 1:N
        T(j, 1) = (T(j, 2) + (h * dx / lambda) * Tinf) / (1 + (h * dx / lambda));
    end
end