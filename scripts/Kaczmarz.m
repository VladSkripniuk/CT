function F = Kaczmarz(A,g,lambda_schedule,num_iter,verbosity, version, tag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = size(A, 1);

F = zeros(size(A, 2),1);

saveplots_after = round(N / 10, 0);

mkdir(sprintf('../pics/%s', tag));

fileID = fopen(sprintf('../pics/%s/log.txt', tag),'w');
fprintf("before step 1 MSE: %.5f\n", mean((A*F-g).^2));
fprintf(fileID, "before step 1 MSE: %.5f\n", mean((A*F-g).^2));

for k=1:num_iter
    
    lambda = lambda_schedule(k);
    permutation_of_row_indices = randperm(N);
    
    for step_index=1:N
        
        if version == "randomized"
            row_index = permutation_of_row_indices(step_index);
        else
            row_index = step_index;
        end
        
        aj = A(row_index,:);

        F = F + lambda * (g(row_index) - dot(aj, F)) * aj' / norm(aj)^2;
        F(F > 1) = 1;
        F(F < 0) = 0;
        
        if strcmp(verbosity, "saveplots") && mod(step_index, saveplots_after) == 0
   
            img = reshape(F, round(sqrt(size(F, 1))), round(sqrt(size(F, 1))));
            
            resultsfigure = figure('visible', 'off');
            imshow(img,[min(F) max(F)]);
            print(resultsfigure, '-dpng', sprintf('../pics/%s/%d.png', tag, k*1000000+step_index), '-r300');
            close;
        end
    end
    fprintf("after step %d MSE: %.5f\n", k, mean((A*F-g).^2));
    fprintf(fileID, "after step %d MSE: %.5f\n", k, mean((A*F-g).^2));
end

fclose(fileID);

end

