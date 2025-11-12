clc;
clear all;

%% Solving the governing equation using line by line solver by making tdma structure

% Parameters 
thermal_diffusivity = 5e-6;
h = 250;
lambda = 5;
T_inf = 573;
T_initial = 300;
length = 0.08;
n = 81;
delta_x = length/(n-1);
delta_t = 0.1;
tolerance = 1e-6;

time_instance = 6;
convergence_time = 1828.70;

%% using implicit scheme to solve the differential equation 

% Initialize temperature field initially
% initial temperature condition
T = T_initial*ones(n,n);

% specifying the boundary conditions 
T(1, 1:n) = 423;     
T(1:n, n) = 473;      
T(n, 1:n) = 323;  

for i = 2:n-1
    T(i,1) = (lambda * T(i,2) + h * delta_x * T_inf) / (lambda + h * delta_x);
end

start_time = tic;
target_times = linspace(0, convergence_time, time_instance);

% store temperature in matrices
vertical_data = zeros(time_instance, n); % 6 rows and 81 columns 1 row for 1 time-temperature
horizontal_data = zeros(time_instance, n);  
captured_times = zeros(1, time_instance); % Initialize captured_times

% Initial condition is stored in 1st row
vertical_data(1, :) = T(:, 41);
horizontal_data(1, :) = T(41, :);
captured_times(1) = 0;
current_index = 2; % maintaining the next row index that needs to be updated

[T, converged_time, vertical_data, horizontal_data, captured_times] = line_by_line_solver(T, n, delta_x, delta_t, h, lambda, T_inf, tolerance, thermal_diffusivity, target_times, vertical_data, horizontal_data, captured_times, current_index, time_instance);

end_time = toc(start_time);

fprintf('Steady state reached at time: %.2f seconds\n', converged_time);
fprintf('Computational time: %.4f seconds\n', end_time);

%% Results Plotting
fprintf('\nPlotting results...\n');

% Create position vectors for plotting
y_positions = (0:n-1) * delta_x * 1000; % Convert to mm
x_positions = (0:n-1) * delta_x * 1000; % Convert to mm

% Plot 1: Vertical centerline

figure;
hold on;
for i = 1:time_instance
    plot(y_positions,vertical_data(i, :), 'LineWidth', 2);
end
xlabel('Temperature (K)');
ylabel('Y Position (mm)');
title('Temperature Variation along Vertical Centerline (x = 40 mm)');
legend({'0.00', '365.74', '731.48', '1097.22', '1462.96', '1828.70'}, 'Location', 'best')
grid on;
hold off;

% Plot 2: Horizontal centerline 
figure;
hold on;
for i = 1:time_instance
    plot(x_positions, horizontal_data(i, :), 'LineWidth', 2);
end
xlabel('X Position (mm)');
ylabel('Temperature (K)');
title('Temperature Variation along Horizontal Centerline');
legend({'0.00', '365.74', '731.48', '1097.22', '1462.96', '1828.70'}, 'Location', 'best')
grid on;
hold off;

%% contour plot

figure;
contourf(flipud(T), 40);
colorbar;
title('Temperature Field - 81 x 81 Using Implicit Scheme');
xlabel('X');
ylabel('Y');
axis equal tight;

%% making the function from line by line solver

function [T, converged_time, vertical_data, horizontal_data, captured_times] = line_by_line_solver(T, n, delta_x, delta_t, h, lambda, T_inf, tolerance, thermal_diffusivity, target_times, vertical_data, horizontal_data, captured_times, current_index, time_instance)
    max_time_steps = 1e+6;
    
    for k = 1:max_time_steps
        T_old = T;
        
        % Row sweep i is constant
        for i = 2:n-1
            n_unknowns = n - 2;  % making the boundaries known
            
            % Mat A and rhs mat b
            A = zeros(n_unknowns, n_unknowns);
            b = zeros(n_unknowns, 1);
            
            for j_idx = 1:n_unknowns
                j = j_idx + 1;  % Actual j coordinate (j=2 to n-1)
                
                % Main diagonal element
                A(j_idx, j_idx) = (1 + (4*thermal_diffusivity*delta_t)/((delta_x)^2));
                
                % left element
                if j_idx > 1
                    A(j_idx, j_idx-1) = -(thermal_diffusivity*delta_t)/((delta_x)^2);
                else
                    % if left element is left boundary
                    A(j_idx, j_idx) = A(j_idx, j_idx) - ((thermal_diffusivity*delta_t)/((delta_x)^2))*(lambda/(lambda + h*delta_x));
                end
                
                % right element
                if j_idx < n_unknowns
                    A(j_idx, j_idx+1) = -(thermal_diffusivity*delta_t)/((delta_x)^2);
                end
                
                % b matrix
                b(j_idx) = T_old(i,j) + ((thermal_diffusivity*delta_t)/((delta_x)^2)) * (T_old(i-1,j) + T_old(i+1,j));
                
                % Boundary contributions to b
                if j_idx == 1
                    % Left boundary condition
                    b(j_idx) = b(j_idx) + ((thermal_diffusivity*delta_t)/((delta_x)^2))*((h*delta_x*T_inf)/(lambda + h*delta_x));
                elseif j_idx == n_unknowns
                    % Right boundary
                    b(j_idx) = b(j_idx) + ((thermal_diffusivity*delta_t)/((delta_x)^2)) * T(i, n);
                end
            end
            
            % Solve using your TDMA function
            if n_unknowns > 0
                T_row = tdma(A, b);
                T(i, 2:n-1) = T_row';
            end
            
            % Update the left boundary with new interior values
            T(i, 1) = (lambda * T(i, 2) + h * delta_x * T_inf) / (lambda + h * delta_x);
        end
        
        % Calculate maximum residual
        max_residual = 0;
        for i = 2:n-1
            for j = 2:n-1
                residual = abs(T(i,j) - T_old(i,j));
                if residual > max_residual
                    max_residual = residual;
                end
            end
        end

        current_time = k * delta_t;

        % Capture data at target times
        if current_index <= time_instance
            if current_time >= target_times(current_index)
                vertical_data(current_index, :) = T(:, 41);
                horizontal_data(current_index, :) = T(41, :);
                captured_times(current_index) = current_time;
                current_index = current_index + 1;
            end
        end
        
        
        % Check for steady state convergence
        if max_residual < tolerance
            converged_time = current_time;
            return;
        end
    end
end

%% tutorial tdma 

function x = tdma(A, b)
    n = length(b);
    x = zeros(n, 1);
    
    % Forward elimination
    for i = 2:n
        factor = A(i, i-1) / A(i-1, i-1);
        A(i, i) = A(i, i) - factor * A(i-1, i);
        b(i) = b(i) - factor * b(i-1);
    end
    
    % Back substitution
    x(n) = b(n) / A(n, n);
    for i = n-1:-1:1
        x(i) = (b(i) - A(i, i+1) * x(i+1)) / A(i, i);
    end
end