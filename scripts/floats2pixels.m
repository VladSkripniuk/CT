function [i, j] = floats2pixels(x, y, res)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
if x > 1 - 1e-6
    x = 1 - 1e-6;
end
if x < -(1 - 1e-6)
    x = -(1 - 1e-6);
end
if y > 1 - 1e-6
    y = 1 - 1e-6;
end
if y < -(1 - 1e-6)
    y = -(1 - 1e-6);
end
i = round(-y*(res+0.5), 0) + res+1;
j = round(x*(res+0.5), 0) + res+1;
end