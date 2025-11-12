clc;
clear all;

%% Validation Parameters

p_exp_psi = [0.9608492, 0.9606546, 0.9600656, 0.9590641, 0.9576204, 0.9556900, ...
0.9532128, 0.9501103, 0.9462820, 0.9416005, 0.9359082, 0.9290103, ...
0.9206690, 0.9105966, 0.8984537, 0.8838449, 0.8663263, 0.8454205, ...
0.8206479, 0.7915765, 0.7578946, 0.7194992, 0.6765875, 0.6297246, ...
0.5798628, 0.5286468, 0.4983042, 0.4686222, 0.4395775, 0.4114938, ...
0.3846425, 0.3592383, 0.3354384, 0.3133435, 0.2930009, 0.2744158, ...
0.2575561, 0.2423654, 0.2287681, 0.2166777, 0.2060031, 0.1966533, ...
0.1885405, 0.1815824, 0.1757052, 0.1708438, 0.1669425, 0.1639555, ...
0.1618473, 0.1605926, 0.1601759];

p_exp_Pa = p_exp_psi * 6894.76; % the data is given in opsi convert to pasacal

% Experimental Mach data
mach_exp = [0.2395428, 0.2401525, 0.2419912, 0.2450882, 0.2494931, 0.2552769, ...
0.2625332, 0.2713805, 0.2819642, 0.2944594, 0.3090735, 0.3260485, ...
0.3456640, 0.3682380, 0.3941266, 0.4237210, 0.4574398, 0.4957166, ...
0.5389794, 0.5876219, 0.6419665, 0.7022194, 0.7684249, 0.8404243, ...
0.9178230, 0.9997698, 1.049278, 1.099540, 1.150450, 1.201643, ...
1.252754, 1.303419, 1.353288, 1.402030, 1.449334, 1.494914, ...
1.538515, 1.579906, 1.618885, 1.655275, 1.688924, 1.719701, ...
1.747494, 1.772210, 1.793771, 1.812114, 1.827187, 1.838950, ...
1.847373, 1.852435, 1.854124];

%% Parameters for isentropic subsonic-supersonic case

n = 51;                    % Number of grid points                
x = linspace(0, 0.254, n); % Axial coordinate
dx = x(2) - x(1);          
gamma = 1.4;              
c = 0.5;                   % CFL number

%% Defining the data as per NPARC geometry data 

a = zeros(1, n); % area 1D matrix formulation
for i = 1:n
    if x(i) < 0.127
        a(i) = 0.0444 - 0.019 * cos((0.2*x(i)/0.0254 - 1) * pi);
    else
        a(i) = 0.0318 - 0.0063 * cos((0.2*x(i)/0.0254 - 1) * pi);
    end
end

% Find throat location (minimum area)
[~, throat] = min(a);
A_throat = a(throat);

%% Boundary conditions for isentropic subsonic-supersonic case
% Low back pressure causes flow to become supersonic

P0 = 6894.76;                % Total pressure at inlet (Pa) not the stangnaton pressure
T0 = 55.56;                  %  Temperature at inlet (K)  
P_back = 1103.16;            % Back pressure (Pa)
fprintf('Inlet total conditions: P0 = %.2f Pa, T0 = %.2f K\n', P0, T0);

%% function to solve the nozzel flow equations 

[sim_time_c, p_ratio,rho_ratio,T_ratio, mach_final, p_final, rho_final, T_final, m_dot_final] = solver(n, x, dx, a, gamma, P0, T0,throat, P_back);

%% Results 

fprintf('Simulation time: %.3g seconds\n', sim_time_c);

fprintf('Mach Number: %.4f\n', mach_final(throat));
fprintf('Pressure Ratio (P/P0): %.4f\n', p_ratio(throat));
fprintf('Density Ratio (ρ/ρ0): %.4f\n', rho_ratio(throat));
fprintf('Temperature Ratio (T/T0): %.4f\n', T_ratio(throat));
fprintf('Average Mass Flow Rate: %.6f kg/s\n', mean(m_dot_final));

figure;
plot(x, a, 'k-', 'LineWidth', 2)
hold on
plot(x(throat), a(throat), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 2)
xlabel('Axial Distance (m)')
ylabel('Area (m²)')
title('Nozzle Geometry - Area')
legend('Area Distribution', 'Throat', 'Location', 'best')
grid on

%% dimensional plotting
figure;
plot(x, p_final, '-r', 'LineWidth', 2)
ylabel('Pressure ( P )')
xlabel('Axial Distance (m)')
grid on
title('Pressure Distribution Along Axis')

figure;
plot(x, rho_final, '-g', 'LineWidth', 2)
ylabel('Density (rho) ')
xlabel('Axial Distance (m)')
grid on
title('Density Distribution ALong Axis')

figure;
plot(x, T_final, '-m', 'LineWidth', 2)
hold on
ylabel('Temperature ( T )')
xlabel('Axial Distance (m)')
grid on
title('Temperature Distribution Along Axis')
%% non dimensional plotting

figure;
plot(x, p_ratio, 'LineWidth', 1.5)
ylabel('Pressure Ratio (P/P_0)')
xlabel('Axial Distance (m)')
grid minor
hold on
legend('Pressure','Location', 'northeast')
title('Non Dimensional Pressure Distribution Along Axial Direction')

figure;
plot(x, T_ratio, 'LineWidth', 1.5)
xlabel('Axial Distance (m)')
ylabel('Temperature Ratio (T/T_0)')
grid on
hold on
legend('Temperature', 'Location', 'northeast')
title('Temperature Distribution')

figure;
plot(x, rho_ratio, '-g', 'LineWidth', 2)
ylabel('Density Ratio (ρ/ρ_0)')
xlabel('Axial Distance (m)')
grid minor
title('Density Distribution')
%% Figure 5: Mass Flow Rate Distribution
figure;
plot(x, m_dot_final, 'b-o', 'LineWidth', 1.5)
xlabel('Axial Distance (m)')
ylabel('Mass Flow Rate (kg/s)')
title('Mass Flow Rate Distribution Along Nozzle')
grid on
hold on
yline(mean(m_dot_final), 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Mean = %.6f kg/s', mean(m_dot_final)))
legend('Mass flow rate', 'Mean', 'Location', 'best')

%% Validation Plots

figure;
plot(x, p_final, 'b-', 'LineWidth', 2);
hold on;
plot(x, p_exp_Pa, 'ro','LineWidth', 2 );
xlabel('Axial Distance (m)');
ylabel(' Pressure (Pa)');
title('Comparison of Pressure Distribution: Numerical vs Experimental');
legend('Numerical Simulation', 'Experimental', 'Location', 'best');
grid on;

figure;
plot(x, mach_final, 'k-', 'LineWidth', 2); 
hold on;
plot(x, mach_exp, 'gs','LineWidth', 2 );
xlabel('Axial Distance (m)');
ylabel('Mach Number');
title('Comparison of Mach Number Distribution: Numerrical vs Experimental');
legend('Numerical Simulation', 'Experimental', 'Location', 'best');
grid on;

%% analytical Solution 
P_P0 = (1 + (gamma - 1)/2 .* mach_final.^2).^(-gamma/(gamma-1));
rho_rho0 = (1 + (gamma - 1)/2 .* mach_final.^2).^(-1/(gamma-1));
T_T0 = (1 + (gamma - 1)/2 .* mach_final.^2).^(-1);

figure;
plot(x, p_ratio, 'b-', 'LineWidth', 2);
hold on;
plot(x, P_P0, 'ro','LineWidth', 2 );
xlabel('Axial Distance (m)');
ylabel(' Pressure Ratio');
title('Comparison of Pressure Distribution: Numerical vs Analytical');
legend('Numerical Simulation', 'Analytical', 'Location', 'best');
grid on;

figure;
plot(x, rho_ratio, 'k-', 'LineWidth', 2); 
hold on;
plot(x, rho_rho0, 'gs','LineWidth', 2 );
xlabel('Axial Distance (m)');
ylabel('Density Ratio');
title('Comparison of ratio of density: Numerical vs Analytical');
legend('Numerical Simulation', 'Analytical', 'Location', 'best');
grid on;

figure;
plot(x, T_ratio, 'k-', 'LineWidth', 2); 
hold on;
plot(x, T_T0, 'gs','LineWidth', 2 );
xlabel('Axial Distance (m)');
ylabel('Temperature Ratio');
title('Comparison of temperature ratio: Numerical vs Analytical');
legend('Numerical Simulation', 'Analytical', 'Location', 'best');
grid on;
%% function to solve the case
function [sim_time_c,p_ratio,rho_ratio,T_ratio, mach_final, p_final, rho_final, T_final, m_dot_final] = solver(n, x, dx, a, gamma, P0, T0,throat, P_back)

    R = 287;
    rho0 = P0 / (R * T0);
    
    % initial conditions by trial and error but doesnt depend for steady state solution 
    rho_c = ones(1,n) * rho0;
    v_c = linspace(50, 200, n);
    T_c = ones(1,n) * T0;

    % Initialize conservative variables
    Q1 = rho_c .* a;
    Q2 = rho_c .* a .* v_c;
    Q3 = rho_c .* a .* ( (T_c*R)/(gamma-1) + 0.5*v_c.^2 );
    
    % Pre-allocate arrays for throat convergence from 1 to 500000 whihc is
    % total numer of itrations 

    % iterations 
    th_mach_c = zeros(1, 50000);
    th_press_c = zeros(1, 50000);
    th_temp_c = zeros(1, 50000);
    th_rho_c = zeros(1, 50000);
    
    % Time marching loop
    tic;
    for k = 1:50000
        
        % Store old values
        Q1_old = Q1;
        Q2_old = Q2;
        Q3_old = Q3;
        
        % Calculate primitive variables for convergence 
        rho_temp = Q1 ./ a;
        v_temp = Q2 ./ Q1;
        e_temp = Q3 ./ Q1;
        T_temp = (gamma-1)/R * (e_temp - 0.5*v_temp.^2);
        p_temp = rho_temp .* R .* T_temp;
        speed_of_sound_local = sqrt(gamma * R * T_c);
        
        % Calculate fluxes
        F1 = Q2;
        F2 = (Q2.^2 ./ Q1) + p_temp .* a;
        F3 = (Q3 + p_temp .* a) .* v_temp;
        
        % Predictor step
        dQ1dt_p = zeros(1, n);
        dQ2dt_p = zeros(1, n);
        dQ3dt_p = zeros(1, n);
        
        for i = 2:n-1
            da_dx = (a(i+1) - a(i-1)) / (2*dx); % centra difference for area
            dlnA_dx = da_dx / a(i); 
     
            last_term = Q3(i) - (Q2(i)^2) / (2 * Q1(i));
            
            J2 = (gamma - 1) * dlnA_dx * last_term;
            
            dQ1dt_p(i) = -(F1(i+1) - F1(i)) / dx;
            dQ2dt_p(i) = -(F2(i+1) - F2(i)) / dx + J2;  % Source term addition
            dQ3dt_p(i) = -(F3(i+1) - F3(i)) / dx;
        end

        % timestep based on cfl from mc cromack's scheme 
        dt = (0.5) * dx/(max(speed_of_sound_local) + max(v_c));
        
        % Predicted values
        Q1_star = Q1 + dQ1dt_p * dt;
        Q2_star = Q2 + dQ2dt_p * dt;
        Q3_star = Q3 + dQ3dt_p * dt;
        
        % Calculate star primitive variables
        rho_star = Q1_star ./ a;
        v_star = Q2_star ./ Q1_star;
        e_star = Q3_star ./ Q1_star;
        T_star = (gamma-1)/R * (e_star - 0.5*v_star.^2);
        p_star = rho_star .* R .* T_star;
        
        % Star fluxes
        F1_star = Q2_star;
        F2_star = (Q2_star.^2 ./ Q1_star) + p_star .* a;
        F3_star = (Q3_star + p_star .* a) .* v_star;
        
        % Corrector step (backward differences)
        dQ1dt_c = zeros(1, n);
        dQ2dt_c = zeros(1, n);
        dQ3dt_c = zeros(1, n);
        
        for i = 2:n-1
            da_dx = (a(i+1) - a(i-1)) / (2*dx);
            dlnA_dx = da_dx / a(i);
            last_star = Q3_star(i) - (Q2_star(i)^2) / (2 * Q1_star(i));
            J2_star = (gamma - 1) * dlnA_dx * last_star;
            
            dQ1dt_c(i) = -(F1_star(i) - F1_star(i-1)) / dx;
            dQ2dt_c(i) = -(F2_star(i) - F2_star(i-1)) / dx + J2_star;
            dQ3dt_c(i) = -(F3_star(i) - F3_star(i-1)) / dx;
        end
        
        % Average derivatives
        dQ1dt_avg = 0.5 * (dQ1dt_p + dQ1dt_c);
        dQ2dt_avg = 0.5 * (dQ2dt_p + dQ2dt_c);
        dQ3dt_avg = 0.5 * (dQ3dt_p + dQ3dt_c);
        
        % Update solution
        Q1 = Q1_old + dQ1dt_avg * dt;
        Q2 = Q2_old + dQ2dt_avg * dt;
        Q3 = Q3_old + dQ3dt_avg * dt;
        
        % BOUNDARY CONDITIONS
        % Calculate current primitive variables
        v_c = Q2 ./ Q1;
        e_c = Q3 ./ Q1;
        T_c = (gamma-1)/R * (e_c - 0.5*v_c.^2);
        
        % INLET
        v_in = 2*v_c(2) - v_c(3);
        T_in = T0 - v_in^2 / (2 * (R*gamma/(gamma-1)));
        p_in = P0 * (T_in/T0)^(gamma/(gamma-1));
        rho_in = p_in / (R * T_in);
        
        % Apply inlet boundary conditions for the conservative varaibels
        Q1(1) = rho_in * a(1);
        Q2(1) = rho_in * a(1) * v_in;
        Q3(1) = rho_in * a(1) * ((T_in*R)/(gamma-1) + 0.5*v_in^2);
        
        % OUTLET extrapolation
        Q1(n) = 2*Q1(n-1) - Q1(n-2);
        Q2(n) = 2*Q2(n-1) - Q2(n-2);
        Q3(n) = 2*Q3(n-1) - Q3(n-2);
        
        % Update primitive variables after boundary conditions
        rho_c = Q1 ./ a;
        v_c = Q2 ./ Q1;
        e_c = Q3 ./ Q1;
        T_c = (gamma-1)/R * (e_c - 0.5*v_c.^2);
        
        % Calculate derived quantities
        p_c = rho_c .* R .* T_c;
        mach_no = v_c ./ speed_of_sound_local;
        m_dot = rho_c .* a .* v_c;
        
        % Non-dimensionalize the primitive variables for results plotting
        p_ratio = p_c / P0;
        rho_ratio = rho_c / rho0;
        T_ratio = T_c / T0;
        
        % Store throat values for convergence history
        th_mach_c(k) = mach_no(throat);
        th_press_c(k) = p_ratio(throat);
        th_rho_c(k) = rho_ratio(throat);
        th_temp_c(k) = T_ratio(throat);
        
        % Convergence check
        if k > 1 % check after 1st 
            mach_change = abs(th_mach_c(k) - th_mach_c(k-1));
            press_change = abs(th_press_c(k) - th_press_c(k-1));
            if mach_change < 1e-8 && press_change < 1e-8
                fprintf('Converged at iteration %d\n', k);
                break;
            end
        end
    end
    
    % Final calculations
    mach_final = mach_no;
    p_final = p_c;
    rho_final = rho_c;
    T_final = T_c;
    m_dot_final = m_dot;
    sim_time_c = toc;
end