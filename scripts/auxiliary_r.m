function [r] = auxiliary_r(s)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
s = s / pi;
r = sinc(s) - 0.5*(sinc(0.5*s).^2);
r(find(abs(s) < 1e-4)) = 0.5;
end

