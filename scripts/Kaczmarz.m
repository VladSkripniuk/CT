function F = Kaczmarz(A,g,lambda,num_iter)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = size(A, 1);

% F = rand(size(A, 2),1)-0.5;
F = zeros(size(A, 2),1)-0.5;
F_0 = F;

for k=1:num_iter
    for j=1:N
        aj = A(j,:);
%         disp(size(aj));
%         disp(size(F));
%         fprintf("offset: %.5f\n", (g(j) - dot(aj, F)));
%         fprintf("dF: %.5f\n", norm(lambda * (g(j) - dot(aj, F)) / (norm(aj)^2) * aj));
        F = F + lambda * (g(j) - dot(aj, F)) / (norm(aj)^2) * aj';
    end
    fprintf("MSE: %.5f\n", norm(A*F-g))

end

fprintf("change: %.5f\n", norm(F-F_0))
fprintf("change AF: %.5f\n", norm(A*F-A*F_0))

end

