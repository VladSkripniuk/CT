function radon_transform = CalculateRadonTransformCircle(circledata, s, omega)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    xcenter = circledata(1);
    ycenter = circledata(2);
    radius = circledata(3);
    a = circledata(4);
    b = circledata(5);
    a = a * radius;
    b = b * radius;
    intensity = circledata(6);
    alpha = circledata(7);

    xM = [xcenter ycenter];
    
    ro = sqrt(dot([a^2 b^2], (rotate_vector(omega, -alpha)).^2));
    
    if abs(s - dot(omega, xM)) <= ro
        radon_transform = intensity * 2*abs(a*b)/(ro^2)*sqrt(ro^2-(s - dot(omega, xM))^2);
    else
        radon_transform = 0;
    end             
end

