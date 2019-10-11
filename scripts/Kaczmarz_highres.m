function F = Kaczmarz_highres(M, p, q, g,lambda_schedule,num_iter,verbosity,version, tag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = size(g, 1);

F = zeros(M,1);

saveplots_after = round(N / 10, 0);

mkdir(sprintf('../pics/%s', tag));

% fileID = fopen(sprintf('../pics/%s/log.txt', tag),'w');
% fprintf("before step 1 MSE: %.5f\n", mean((A*F-g).^2));
% fprintf(fileID, "before step 1 MSE: %.5f\n", mean((A*F-g).^2));

for k=1:num_iter
    
    lambda = lambda_schedule(k);
    permutation_of_row_indices = randperm(N);
    
    for step_index=1:N
        
        if version == "randomized"
            row_index = permutation_of_row_indices(step_index);
        else
            row_index = step_index;
        end

        j = mod(row_index-1, 2*q+1) - q;
        i = floor((row_index-1) / (2*q+1)) + 1;

        phi = (i-1)*pi/p;
        s = j / q;
        
        A_row = CalculateApproximateIntersections(phi, s, M);
        A_row = reshape(A_row, 1, M);

        intersections_with_boundary = CalculatePixelIntersections(phi, s, 1);
%         index = 2;
%         while index < size(intersections_with_boundary, 1) ...
%                 && norm(intersections_with_boundary(1, :)-intersections_with_boundary(index, :)) < 1e-2
%             index = index + 1;
%         end
%         
        L = norm(intersections_with_boundary(1, :)-intersections_with_boundary(size(intersections_with_boundary,1), :));
        L = L / (2/sqrt(M));
        
%         fprintf("i: %d, j: %d, row_index: %d\n", i, j, row_index);
%         fprintf("phi: %.5f, s: %.5f, L: %.5f, A_sum:%.5f\n", phi, s, L, sum(A_row));
%         fprintf("g: %.5f, AF: %.5f, g/L: %.5f, AF/sumA:%.5f\n", g(row_index), dot(A_row, F), g(row_index)/L, dot(A_row, F) / sum(A_row));
%         fprintf("update: %.5f\n", (g(row_index)/L - dot(A_row, F) / sum(A_row)));
        F = F + lambda * (g(row_index)/L - dot(A_row, F) / sum(A_row)) * A_row';
%         F = F + lambda * (g(row_index) - dot(A_row, F)) / sum(A_row) * A_row';

%         fprintf("g: %.5f, AF: %.5f, g/L: %.5f, AF/sumA:%.5f\n", g(row_index), dot(A_row, F), g(row_index)/L, dot(A_row, F) / sum(A_row));

        if max(max(F)) > 1
            fprintf("i: %d, j: %d, row_index: %d\n", i, j, row_index);
            fprintf("phi: %.5f, s: %.5f, L: %.5f, A_sum:%.5f\n", phi, s, L, sum(A_row));
            fprintf("g: %.5f, AF: %.5f, g/L: %.5f, AF/sumA:%.5f\n", g(row_index), dot(A_row, F), g(row_index)/L, dot(A_row, F) / sum(A_row));
            fprintf("update: %.5f\n", (g(row_index)/L - dot(A_row, F) / sum(A_row)));
            lkjlfakjd
        end

        F(F > 1) = 1;
        F(F < 0) = 0;
       
%         disp(max(max(F)));
        
        if strcmp(verbosity, "saveplots") && mod(step_index, saveplots_after) == 0
   
            img = reshape(F, round(sqrt(size(F, 1))), round(sqrt(size(F, 1))));
            
            resultsfigure = figure('visible', 'off');
            imshow(img,[min(F) max(F)]);
            print(resultsfigure, '-dpng', sprintf('../pics/%s/%d.png', tag, k*1000000+step_index), '-r300');
            close;
        end
    end
    
%     fprintf("after step %d MSE: %.5f\n", k, mean((A*F-g).^2));
%     fprintf(fileID, "after step %d MSE: %.5f\n", k, mean((A*F-g).^2));
end

% fclose(fileID);

end

