function sinogram = GenerateMeasuredData(p, q, objects)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
nobjects = size(objects, 1);
sinogram = zeros(p, 2*q+1);
delta_phi = pi/p;
delta_s = 1/q;

for obj_index=1:nobjects
    for i = 1:p
        phi = (i-1)*delta_phi;
        omega = [cos(phi) sin(phi)];
        for j = -q:q
            s = delta_s*j;
            sinogram(i, j+q+1) = sinogram(i, j+q+1) + CalculateRadonTransformCircle(objects(obj_index,:), s, omega);
        end
    end
end
end

