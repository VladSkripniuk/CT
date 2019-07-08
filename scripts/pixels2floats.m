function [x, y] = pixels2floats(i, j, res)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
x = -1 + (j-1) * 2 / (2*res+1) + 1 / (2*res+1);
y = 1 - (i-1) * 2 / (2*res+1) - 1 / (2*res+1);
end

