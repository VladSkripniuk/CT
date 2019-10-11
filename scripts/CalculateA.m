function A = CalculateA(p,q,M,angle_restr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('angle_restr','var')
    % third parameter does not exist, so default it to something
    angle_restr = 0;
end

A = zeros(p*(2*q+1), M);
filename = sprintf("A_p%d_q%d_M%d_angle%d.mat", p, q, M, round(1000*angle_restr, 0));
path = sprintf("../matrices/%s", filename);

max_angle = pi - angle_restr;
min_angle = angle_restr;
delta_phi = (max_angle-min_angle) / p;

if isfile(path)
    mat = load(path, 'A');
    A = mat.A;
else
    for i=1:p
        for j=-q:q
%             phi = (i-1)*pi/p;
            phi = min_angle + (i-1)*delta_phi;
            s = j / q;

            pixel_intersections = CalculatePixelIntersections(phi, s, M);

            A((i-1)*(2*q+1)+j+q+1,:) = CalculateDistances(pixel_intersections, M);
        end
    end

    % save as .mat file
    save(path, 'A', '-v7.3');
end

end