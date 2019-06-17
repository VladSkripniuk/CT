function phantom = PhantomCircle(parameter,res)
%PhatomCircle Summary of this function goes here
%   Detailed explanation goes here
phantom = zeros(2*res+1,2*res+1);

nobjects = size(parameter, 1);
for obj_index=1:nobjects
    xcenter = parameter(obj_index, 1);
    ycenter = parameter(obj_index, 2);
    radius = parameter(obj_index, 3);
    a = parameter(obj_index, 4);
    b = parameter(obj_index, 5);
    intensity = parameter(obj_index, 6);
    alpha = parameter(obj_index, 7);

    %     https://math.stackexchange.com/questions/91132/how-to-get-the-limits-of-rotated-ellipse
    xleft = xcenter - sqrt(a^2 * radius^2 * cos(alpha)^2 + b^2 * radius^2 * sin(alpha)^2);
    xright = xcenter + sqrt(a^2 * radius^2 * cos(alpha)^2 + b^2 * radius^2 * sin(alpha)^2);
    ytop = ycenter + sqrt(a^2 * radius^2 * sin(alpha)^2 + b^2 * radius^2 * cos(alpha)^2);
    ybottom = ycenter - sqrt(a^2 * radius^2 * sin(alpha)^2 + b^2 * radius^2 * cos(alpha)^2);

    [itop, jtop] = floats2pixels(xleft, ytop, res);
    [ibottom, jbottom] = floats2pixels(xright, ybottom, res);

    for i = itop:ibottom
        for j = jtop:jbottom
            
            [x,y] = pixels2floats(i,j,res);
            x = x - xcenter;
            y = y - ycenter;
            if ((x*cos(alpha)-y*sin(alpha))/a)^2+((x*sin(alpha)+y*cos(alpha))/b)^2 <= radius^2
                phantom(i,j) = phantom(i,j) + intensity;
            end
        end
    end
end
%save('phantom.mat', phantom);

end