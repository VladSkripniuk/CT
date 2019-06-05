function [i, j] = floats2pixels(x, y, res)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
i = round(-y*res, 0) + res+1;
j = round(x*res, 0) + res+1;
end