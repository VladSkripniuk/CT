function [v] = CalculateConvolution(g, p, q, wb)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
v = zeros(size(g));
for i = 1:p
    for k = -q:q
        wb_shifted = zeros(size(wb));
        %left_bound = min(1, 1+k);
        %right_bound = max(2*q+1, 2*q+1+k);
        %wb_shifted(left_bound:right_bound)= wb(left_bound:right_bound);
        wb_shifted(max(1,1+k):min(2*q+1,2*q+1+k)) = wb(max(1,1-k):min(2*q+1,2*q+1-k));
        v(i, q+1+k) = sum(wb_shifted .* g(i,:)) / q;
    end
end
end

