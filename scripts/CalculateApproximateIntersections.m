function A_row = CalculateApproximateIntersections(phi, s, M)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm

eps=1e-4;

res = (sqrt(M)-1) / 2;
delta = 2 / (2*res+1);
A_row = zeros(2*res+1);

k = cot(phi);
omega = [cos(phi) sin(phi)];

% disp(omega);
% disp(k);

transposed_flag = 0;

if abs(k) > 1
    transposed_flag = 1;
    k = 1/k;
    omega = [omega(2) omega(1)];
end


% disp(omega);
% disp(k);

for j=1:(2*res+1)
    x = -1 + delta/2 + (j - 1) * delta;
    y = (s - omega(1)*x)/omega(2);
    if abs(y) > 1+ eps
        continue;
    end
    
    
%     disp(x);
%     disp(y);
    
    i = round(-y*res, 0) + res+1;
    
%     disp(i);
%     disp(j);
    
    A_row(i,j) = 1;
end

if transposed_flag == 1
    A_row = rot90(A_row, 2)';
end

end