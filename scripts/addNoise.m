function A = addNoise(A, noiselevel)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
delta = 100; % in percent
N = 1+0.02*delta*(rand(size(A))-0.5);
A = A.*N;
end

