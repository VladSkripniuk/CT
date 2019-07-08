function F = Kaczmarz(A,g,lambda,num_iter,verbosity)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = size(A, 1);

F = zeros(size(A, 2),1);

% labmda0 = lambda;

for k=1:num_iter
    
    for i=1:N
        
% %         lambda = labmda0 / sqrt(((k-1)*N+i)/500);
%         curve((k-1)*N+i) = norm(A*F - g);
        j = randi(N, 1);
        aj = A(j,:);
%         disp(size(aj));
%         disp(size(F));
%         fprintf("offset: %.5f\n", (g(j) - dot(aj, F)));
%         fprintf("dF: %.5f\n", norm(lambda * (g(j) - dot(aj, F)) / (norm(aj)^2) * aj));
%         if norm(aj) < 1e-3
%             continue;
%         end
%       
%         if norm(lambda * (g(j) - dot(aj, F)) * aj' / norm(aj)^2) > 3
%             fprintf("iter: %d, update: %.5f\n", i, norm(lambda * (g(j) - dot(aj, F)) * aj' / norm(aj)^2));
%             fprintf("g: %.5f, diff: %.5f, a: %.5f\n", g(j), g(j) - dot(aj, F), norm(aj));
%         end
%         F0 =F;
        F = F + lambda * (g(j) - dot(aj, F)) * aj' / norm(aj)^2;
        F(F > 1) = 1;
        F(F < 0) = 0;
        if strcmp(verbosity, "saveplots") && mod(i, 1010) == 0
   
            img = reshape(F, round(sqrt(size(F, 1))), round(sqrt(size(F, 1))));
            
            resultsfigure = figure;
            imshow(img,[min(F) max(F)]);
            print(resultsfigure, '-dpng', sprintf('../pics/kaczmarz/%d.png', k*1000000+i), '-r300');
            close;
        end
    end
    fprintf("MSE: %.5f\n", mean((A*F-g).^2))
    lambda = lambda / 2;
end

% fprintf("change: %.5f\n", norm(F-F_0))
% fprintf("change AF: %.5f\n", norm(A*F-A*F_0))

end

