%% Make the gaussian function to combine everything

function gauss = gaussian_elimination(A,b)

size_of_matrix = size(A);
row = size_of_matrix(1);
col = size_of_matrix(2);

for j = 1:col-1 % choose 1st to col-1 columnms for performing the elimination
    for i = j+1:row % number of rows we need to perform 
        factor = A(i,j)/A(j,j);
        for k = j:row % then update for the row
            A(i,k) = A(i,k) - factor*A(j,k);
        end
        b(i) = b(i) - factor*b(j);
    end
end

% performing the backward substitution
solution_mat_final = zeros(row,1);
solution_mat_final(col) = b(col)/A(row,col);

for i = row-1:-1:1
    subtracting_amount = 0;
    for j = col:-1:i+1
        subtracting_amount = subtracting_amount + (solution_mat_final(j)*A(i,j));
    end
    solution_mat_final(i) = (b(i) - subtracting_amount)/A(i,i);
end

gauss = solution_mat_final;

end

