function [x, y] = pixels2floats(i, j, res)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
x = (j-1) / res -1;
y = -((i-1) / res -1);
end

