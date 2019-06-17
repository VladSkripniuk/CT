function [wb] = CalculateFilter(s, b, filterchoice)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if strcmp(filterchoice, 'Ram-Lak')
    wb = b^2 /(4 * pi ^2) * auxiliary_r(b * s);
elseif strcmp(filterchoice, 'Cosine')
    wb = b^2 / (8* pi^2)*(auxiliary_r(pi/2-b*s) + auxiliary_r(pi/2+b*s)); 
elseif strcmp(filterchoice, 'Shepp-Logan')
    wb = b^2 /(2 * pi ^3) * r_hat(b * s);
elseif strcmp(filterchoice, 'no filter')
    wb = zeros(size(s));
    wb(find(abs(s) < 1e-4)) = 1;
end
end

