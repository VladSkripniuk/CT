function [r_hat] = r_hat(s)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
r_hat = (pi / 2 -s.*sin(s))./((pi/2)^2-s.^2);
r_hat(find(abs(abs(s)-pi/2) < 1e-4)) = 1/pi;
end

