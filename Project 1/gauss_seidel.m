function [seidel , iteration] = gauss_seidel(T_left,T_right, tolerance, n)
    x = zeros(n,1); % initial guess vector
    x(1) = T_left;
    x(n) = T_right;

    residue = 1000;
    iteration = 0;
    
    while residue > tolerance
        iteration = iteration + 1;
        xold = x;
        for i = 2:n-1
            x(i) = (x(i+1) + x(i-1)) / 2;
        end
        residue = max(abs(xold - x));
    end

    seidel = x;
end