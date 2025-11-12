clc;
clear all;

%% Solving the governing equation using ADI approach 

% Initializing the boundary conditions
h = 250;                % W/m^2K
lambda = 5;             % W/m-k
Tinf = 573;               %K
L = 0.08;                % length of domain in m 

x_grid = 21; % grid points in x direction
y_grid = 41;     % grid points in y direction
max_iter = 100000;
tolerance = 1e-6;

dx = L/(x_grid-1); % grid space 
dy = L/(y_grid-1);
    
T = zeros(x_grid, y_grid); %  initial guess

% Apply boundary conditions

T(:, 1) = 323; % Bottom
T(:, end) = 423; % Top
T(end, :) = 473; % Right

% Initial guess for left boundary (convective/Neumann BC)
for j = 1:y_grid
    T(1, j) = (lambda * T(2, j) + h * dx * Tinf) / (lambda + h * dx);
end

% Start cpu time
tic;

[T, residuals, iter] = adi_method(T, x_grid, y_grid, dx, dy, h, lambda, Tinf, max_iter, tolerance);

% Recorded time
cpu_time = toc;
fprintf('Grid 21 x 41: %d iterations, %.4f seconds\n', iter, cpu_time);

%% Plot temperature contour for 21x41 grid

figure;
[X, Y] = meshgrid(linspace(0, L, x_grid), linspace(0, L, y_grid));
T_plot = T'; % this is done becuase the temperature field was coming inverted

contourf(X, Y, T_plot, 40);
colorbar;
title('Temperature Field - 21x41 Grid (ADI Method)');
xlabel('x (m)');
ylabel('y (m)');
axis equal;

%% Visualize the convergence history

figure;
plot(1:length(residuals), residuals, 'LineWidth', 1.5);
xlabel('Number of Iterations');
ylabel('Residual');
title('Convergence History for 21x41 Grid (ADI Method)');
grid on;
