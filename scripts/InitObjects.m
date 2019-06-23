% Initialize the Objects
function parameters = InitObjects()
%PhatomCircle Summary of this function goes here
%   Detailed explanation goes here
% x, y, r, a, b
% parameters = [0 0 0.5 0.5 2 128 pi/4;
%                0.7 -0.7 0.25 1 1 255 0.0;
%                0 0 0.5 2 0.5 128 0.0];

parameters = [0     0       1   0.69    0.92    1       0;
              0     -0.0184 1   0.6624  0.874   -0.8    0;
              0.22  0       1   0.11    0.31    -0.2    -18/180*pi;
              -0.22 0       1   0.16    0.41    -0.2    18/180*pi;
              0     0.35    1   0.21    0.25    0.1     0;
              0     0.1     1   0.046   0.046   0.1     0;
              0     -0.1    1   0.046   0.046   0.1     0;
              -0.08 -0.605  1   0.046   0.023   0.1     0;
              0     -0.605  1   0.023   0.023   0.1     0;
              0.06  -0.605  1   0.023   0.046   0.1     0;
              ];
          
end

