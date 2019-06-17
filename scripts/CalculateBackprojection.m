function fFBI = CalculateBackprojection(convolution, p, q, res)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
delta_phi = pi/p;
h = 1/q;
fFBI = zeros(2*res+1,2*res+1); 


for x = -1:(1/res):1
    for y = -1:(1/res):1
        f = 0;
        for i = 1:p
            phi = (i-1)*delta_phi;
            omega = [cos(phi) sin(phi)];
            s = dot(omega, [x y]);
            disp('s');
            disp(s);
            
            k = floor(s/h);
            u = s/h - k;
            
            f = f + 2*pi/p * ((1-u) * convolution(i,k+q+1) + u*convolution(i,k+1+q+1));
        end
            
        x_p = floats2pixels(x);
        y_p = floats2pixels(y);
        
        fFBI(x_p, y_p) = f;
        
        
    end
end    

end