function [T, residuals, iter] = line_by_line_solver(T, N, dx, h, lambda, Tinf, max_iter, tolerance)
    residuals = zeros(max_iter, 1);
    
    for iter = 1:max_iter
        T_old = T;
        max_residual = 0;
        
        % Row sweep i is constant
        for i = 2:N-1
            n_unknowns = N - 2;  % making the boundries known
            
            % Mat A and rhs mat b
            A = zeros(n_unknowns, n_unknowns);
            b = zeros(n_unknowns, 1);
            
            for j_idx = 1:n_unknowns
                j = j_idx + 1;  % Actual j coordinate (j=2 to N-1)
                
                % Main diagonal ekement
                A(j_idx, j_idx) = -2 * ((lambda / (dx^2)) + (lambda / (dx^2)));
                
                % left element
                if j_idx > 1
                    A(j_idx, j_idx-1) = lambda / (dx^2);
                else
                    % if left element is left boundary
                    A(j_idx, j_idx) = A(j_idx, j_idx) + (lambda / (dx^2)) * ((lambda/dx) / (lambda/dx + h));
                end
                
                % right element
                if j_idx < n_unknowns
                    A(j_idx, j_idx+1) = (lambda / (dx^2));
                end
                
                % b matrix
                b(j_idx) = -(lambda / (dx^2)) * (T(i-1, j) + T(i+1, j));
                
                % Boundary contributions to b
                if j_idx == 1
                    % Left boundary condition
                    b(j_idx) = b(j_idx) - (lambda / (dx^2)) * (h * Tinf) / (lambda/dx + h);
                elseif j_idx == n_unknowns
                    % Right boundary
                    b(j_idx) = b(j_idx) - (lambda / (dx^2)) * T(i, N);
                end
            end
            
            % Solve using your TDMA function
            if n_unknowns > 0
                T_row = tdma(A, b);
                T(i, 2:N-1) = T_row';
            end
            
            % Update the old entries with new
            T(i, 1) = (lambda/dx * T(i, 2) + h * Tinf) / (lambda/dx + h);
        end
        
        % tolerance checking 
        for i = 2:N-1
            for j = 2:N-1
                residual = abs(T(i,j) - T_old(i,j));
                if residual > max_residual
                    max_residual = residual;
                end
            end
        end
        
        residuals(iter) = max_residual;
        
        if max_residual < tolerance
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