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

    [itop, jtop] = floats2pixels(xcenter - radius*a, ycenter - radius*b, res);
    [ibottom, jbottom] = floats2pixels(xcenter+radius*a, ycenter+radius*b, res);

    for i = itop:ibottom
        for j = jtop:jbottom
            
            [x,y] = pixels2floats(i,j,res);
            if ((x-xcenter)/a)^2+((y-ycenter)/b)^2 <= radius^2
                phantom(i,j) = phantom(i,j) + intensity;
            end
        end
    end
end

end

