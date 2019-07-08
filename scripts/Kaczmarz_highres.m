function F = Kaczmarz_highres(M, p, q, g,lambda,num_iter,verbosity)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = size(g, 1);

F = zeros(M,1);

for k=1:num_iter
    
    permutation_of_row_indices = randperm(N);
    
    for step_index=1:N
        
        row_index = permutation_of_row_indices(step_index);

        j = mod(row_index-1, 2*q+1) - q;
        i = floor(row_index / (2*q+1)) + 1;
%         (i-1)*(2*q+1)+j+q+1

        phi = (i-1)*pi/p;
        s = j / q;
        
        A_row = CalculateApproximateIntersections(phi, s, M);
        A_row = reshape(A_row, 1, M);

        intersections_with_boundary = CalculatePixelIntersections(phi, s, 1);
        index = 2;
        while index < size(intersections_with_boundary, 1) ...
                && norm(intersections_with_boundary(1, :)-intersections_with_boundary(index, :)) < 1e-4
            index = index + 1;
        end
        
%         try
        L = norm(intersections_with_boundary(1, :)-intersections_with_boundary(index, :));
%         catch
%             fprintf("j: %d, i: %d\n", j, i);
%             fprintf("phi: %.5f, s: %.5f\n", phi, s);
%             disp(intersections_with_boundary);
%         end
        L = L / (2/sqrt(M));
        
%         fprintf("i: %d, j: %d, row_index: %d\n", i, j, row_index);
%         fprintf("phi: %.5f, s: %.5f, L: %.5f, A_sum:%.5f\n", phi, s, L, sum(A_row));
%         fprintf("g: %.5f, AF: %.5f, g/L: %.5f, AF/sumA:%.5f\n", g(row_index), dot(A_row, F), g(row_index)/L, dot(A_row, F) / sum(A_row));
%         fprintf("update: %.5f\n", (g(row_index)/L - dot(A_row, F) / sum(A_row)));
        F = F + lambda * (g(row_index)/L - dot(A_row, F) / sum(A_row)) * A_row';

%         fprintf("g: %.5f, AF: %.5f, g/L: %.5f, AF/sumA:%.5f\n", g(row_index), dot(A_row, F), g(row_index)/L, dot(A_row, F) / sum(A_row));
        
        F(F > 1) = 1;
        F(F < 0) = 0;
        
        if strcmp(verbosity, "saveplots") && mod(step_index, 1010) == 0
   
            img = reshape(F, round(sqrt(size(F, 1))), round(sqrt(size(F, 1))));
            
            resultsfigure = figure;
            imshow(img,[min(F) max(F)]);
            print(resultsfigure, '-dpng', sprintf('../pics/kaczmarz/%d.png', k*1000000+step_index), '-r300');
            close;
        end
    end
%     fprintf("MSE: %.5f\n", mean((A*F-g).^2))
    lambda = lambda / 2;
end

end

