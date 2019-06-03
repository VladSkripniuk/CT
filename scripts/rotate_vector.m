function rot_vec = rotate_vector(vec,alpha)
%rotate_vector Rotate vector by a given angle alpha in degrees
vec = vec';
rot = [cos(alpha), -sin(alpha); 
    sin(alpha), cos(alpha)];
rot_vec = rot*vec;
end

