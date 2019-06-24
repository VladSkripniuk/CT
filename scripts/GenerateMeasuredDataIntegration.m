function sinogram = GenerateMeasuredDataIntegration(p, q, phantom)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sinogram = zeros(p, 2*q+1);
delta_phi = pi/p;

for i = 1:p
    for j = -q:q
        phi = (i-1)*delta_phi;
        omega = [cos(phi) sin(phi)];
        s = j/q;
        
        omega_orth = [sin(phi) -cos(phi)];
        
        f = 0;
        for k = -q:q
            pos = s * omega + k/q * omega_orth;
            
            if pos(1)^2 + pos(2)^2 > 1
                continue
            end
            
            [x_pix, y_pix] = floats2pixels(pos(1), pos(2), round((size(phantom, 1)-1)/2));
            f = f + phantom(x_pix, y_pix) / q;
        end
        
        sinogram(i, j+q+1) = f;
        
    end
end

end

