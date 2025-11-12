 clc;
clear all;

%% parameters
thermal_diffusivity = 5e-6;
h = 250;
lambda = 5;
T_inf = 573;
T_right = 473;
T_top = 423;
T_bottom = 323;
T_initial = 300;
length = 0.08;
n = 81;
delta_x = length/(n-1);
delta_t = 0.05;
tolerance = 1e-6;

time_instance = 6;
convergence_time = 915.15;
iterations = 18315;

%% using explicit scheme to solve the differential equations

critical_timestep = (0.5*(delta_x^2))/thermal_diffusivity;
disp(critical_timestep);

% initial temperature condition
T = T_initial*ones(n,n);

% specifying the boundary conditions 
T(1, 1:n) = 423;     
T(1:n, n) = 473;      
T(n, 1:n) = 323;  

for i = 2:n-1
    T(i,1) = (lambda * T(i,2) + h * delta_x * T_inf) / (lambda + h * delta_x);
end

% dividing into 6 equal times
target_times = linspace(0, convergence_time, time_instance);


% store temperature in matrices
vertical_data = zeros(time_instance, n); % 6 rows and 81 columns 1 row for 1 time-temperature
horizontal_data = zeros(time_instance, n);  

% Initial condition is stored in 1st drow
vertical_data(1, :) = T(:, 41);
horizontal_data(1, :) = T(41, :);
captured_times(1) = 0;
current_index = 2; % maintaining the next row index that needs to be updated

starttime = tic;

% we know it wont exceed the iterations but still to ensure 
for k = 1:iterations*iterations
    T_old = T;

    for i = 2:n-1 % for the rows 
        for j = 2:n-1 % for the columns
            T(i,j) = T_old(i,j) + (thermal_diffusivity*delta_t)*((T_old(i+1,j) - 2*T_old(i,j) + T_old(i-1,j))/(delta_x*delta_x) + (T_old(i,j+1) - 2*T_old(i,j) + T_old(i,j-1))/(delta_x*delta_x));
        end
    end

    % using explicit scheme for the left boundary condition as well
    for i = 2:n-1
         T(i,1) = (lambda * T_old(i,2) + h * delta_x * T_inf) / (lambda + h * delta_x);
    end

    current_time = k * delta_t;
    
    % Checking for the particular time 
    if current_index <= time_instance
        if current_time >= target_times(current_index)
            vertical_data(current_index, :) = T(:, 41);
            horizontal_data(current_index, :) = T(41, :);
            captured_times(current_index) = current_time;
            current_index = current_index + 1;
        end
    end
    
    % Check for convergence
    residual = max(max(abs(T - T_old)));
    if residual < tolerance
        fprintf('Converged at = %.2f seconds\n', current_time);
        break;
    end
end

%% Plotting the results

endtime = toc(starttime);
disp(endtime);
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
xlabel('Y Position (mm)');
ylabel('Temperature (K)');
title('Temperature Variation along Vertical Centerline (x = 40 mm)');
legend({'0', '183.03', '366.06', '549.09', '732.12', '915.15'}, 'Location', 'best')
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
legend({'0', '183.03', '366.06', '549.09', '732.12', '915.15'}, 'Location', 'best')
grid on;
hold off;

%% contour plot

figure;
contourf(flipud(T), 40);
colorbar;
title('Temperature Field - 81 x 81 Using Explicit Scheme');
xlabel('X');
ylabel('Y');
axis equal tight;
