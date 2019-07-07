function A = CalculateA(p,q,M)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

A = zeros(p*(2*q+1), M);
filename = sprintf("A_p%d_q%d_M%d.mat", p, q, M);
path = sprintf("../matrices/%s", filename);

if isfile(path)
    mat = load(path, 'A');
    A = mat.A;
else
    for i=1:p
        for j=-q:q
            phi = (i-1)*pi/p;
            s = j / q;

            pixel_intersections = CalculatePixelIntersections(phi, s, M);

            A(i+(j+q)*p,:) = CalculateDistances(pixel_intersections, M);
        end
    end

    % save as .mat file
    save(path, 'A', '-v7.3');
end

end

