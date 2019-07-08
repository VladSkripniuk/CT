function sinogram = GenerateMeasuredData(p, q, objects, angle_restr)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('angle_restr','var')
    % third parameter does not exist, so default it to something
    angle_restr = 0;
end

nobjects = size(objects, 1);
sinogram = zeros(p, 2*q+1);
delta_s = 1/q;

max_angle = pi - angle_restr;
min_angle = angle_restr;
delta_phi = (max_angle-min_angle) / p;

for obj_index=1:nobjects
    for i = 1:p
        phi = min_angle + (i-1)*delta_phi;
        omega = [cos(phi) sin(phi)];
        for j = -q:q
            s = delta_s*j;
            sinogram(i, j+q+1) = sinogram(i, j+q+1) + CalculateRadonTransformCircle(objects(obj_index,:), s, omega);
        end
    end
end
end

