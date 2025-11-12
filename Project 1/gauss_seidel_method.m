function [T, residuals, iter] = gauss_seidel_method(T, N, dx, h, lambda, Tinf, max_iter, tolerance)
    
residuals = zeros(max_iter, 1);

for iter = 1:max_iter

    max_residual = 0;
    
    % Updating interior points 
    for i = 2:N-1
        for j = 2:N-1
            T_new = 0.25 * (T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1));
            residual = abs(T_new - T(i,j));
            if residual > max_residual
                max_residual = residual;
            end
            T(i,j) = T_new;  % old value is updated
        end
    end
    
    % Update left boundary using Gauss-Seidel
    for i = 2:N-1
        T_new = (2*T(i,2) + T(i+1,1) + T(i-1,1) + (2*h*dx*Tinf)/lambda) / (4 + (2*h*dx)/lambda);
        residual = abs(T_new - T(i,1));
        if residual > max_residual
            max_residual = residual;
        end
        T(i,1) = T_new;
    end

    residuals(iter) = max_residual;
    
    % Check convergence
    if max_residual < tolerance
        residuals = residuals(1:iter);
        break;
    end
end

if iter == max_iter
    warning('Maximum iterations reached in the calculations ');
end
end