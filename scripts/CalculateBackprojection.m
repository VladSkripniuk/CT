function fFBI = CalculateBackprojection(convolution, p, q, res)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
delta_phi = pi/p;
h = 1/q;
fFBI = zeros(2*res+1,2*res+1); 


for x = -1:(1/res):1
    for y = -1:(1/res):1
        f = 0;
        
        if x^2 + y^2 > 1
            continue
        end
        
        for i = 1:p
            phi = (i-1)*delta_phi;
            omega = [cos(phi) sin(phi)];
            s = dot(omega, [x y]);
            
            k = floor(s/h);
            u = s/h - k;
            
            if k < q
                f = f + 2*pi/p * ((1-u) * convolution(i,k+q+1) + u*convolution(i,k+1+q+1));
            else
                f = f + 2*pi/p * (1-u) * convolution(i,k+q+1);
            end
        end
            
        [x_p, y_p] = floats2pixels(x, y, res);
        
        fFBI(x_p, y_p) = f;
        
        
    end
end    

end