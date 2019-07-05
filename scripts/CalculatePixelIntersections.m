function intersections = CalculatePixelIntersections(phi, s, M)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

omega = [cos(phi) sin(phi)];
res = (sqrt(M)-1) / 2;

number_of_intersections_with_horizontal = 0;
number_of_intersections_with_vertical = 0;
pixel_intersections_with_horizontal = zeros(2*res+2, 2);
pixel_intersections_with_vertical = zeros(2*res+2, 2);

for y=1:-2/(2*res+1):-1
   %fprintf('y:%.2f', y)
   if (abs(omega(1)) < 1e-4)
       continue
   end
   x = (s - y*omega(2)) / omega(1);
   %fprintf('x:%.2f', x)
   if (abs(x) > 1.004)
       continue
   end
   number_of_intersections_with_horizontal = number_of_intersections_with_horizontal + 1;
   pixel_intersections_with_horizontal(number_of_intersections_with_horizontal, 1) = x;
   pixel_intersections_with_horizontal(number_of_intersections_with_horizontal, 2) = y;
end

for x=1:-2/(2*res+1):-1
   %fprintf('x:%f\n', x)
   if (abs(omega(2)) < 1e-4)
       continue
   end
   y = (s - x*omega(1)) / omega(2);
   %fprintf('y:%f\n', y)
   if (abs(y) > 1.0004)
       continue
   end
   number_of_intersections_with_vertical = number_of_intersections_with_vertical + 1;
   pixel_intersections_with_vertical(number_of_intersections_with_vertical, 1) = x;
   pixel_intersections_with_vertical(number_of_intersections_with_vertical, 2) = y;
end

pixel_intersections_with_horizontal = pixel_intersections_with_horizontal(1:number_of_intersections_with_horizontal,:);
pixel_intersections_with_vertical = pixel_intersections_with_vertical(1:number_of_intersections_with_vertical,:);


%disp(pixel_intersections_with_horizontal)
%disp(pixel_intersections_with_vertical)

if size(pixel_intersections_with_horizontal,1) > 1 && pixel_intersections_with_horizontal(1,1) < pixel_intersections_with_horizontal(2,1)
    pixel_intersections_with_horizontal = flip(pixel_intersections_with_horizontal);
end

intersections = zeros(number_of_intersections_with_horizontal + number_of_intersections_with_vertical, 2);
n = 1;
k = 1;
counter = 0;

for l=1:size(intersections, 1)
    %disp(intersections)
    %fprintf('n: %d', n)
    %fprintf('k: %d', k)
    if n <= size(pixel_intersections_with_horizontal, 1) ...
       && k <= size(pixel_intersections_with_vertical, 1) ...
       && norm(pixel_intersections_with_horizontal(n,:) - pixel_intersections_with_vertical(k,:)) < 1e-2
        n = n+1;
        continue;
    end
    counter = counter + 1;
    
    if n > size(pixel_intersections_with_horizontal, 1)
        intersections(counter,:) = pixel_intersections_with_vertical(k,:);
        k = k+1;
        continue
    end
    if k > size(pixel_intersections_with_vertical, 1)
        intersections(counter,:) = pixel_intersections_with_horizontal(n,:);
        n = n+1;
        continue
    end
    if pixel_intersections_with_horizontal(n,1) > pixel_intersections_with_vertical(k,1)
        intersections(counter,:) = pixel_intersections_with_horizontal(n,:);
        n = n+1;
    else
        intersections(counter,:) = pixel_intersections_with_vertical(k,:);
        k = k+1;
    end
end

intersections = intersections(1:counter,:);

% for j=1:2*res+1
%    x = 1 - (i-1) / res;
%    if (abs(omega(2)) < 1e-4)
%        continue
%    end
%    y = (s - x*omega(1)) / omega(2);
%    if (abs(y) > 1)
%        continue
%    end
%    
%    i = round((x + 1)*res)+1;
%    
%    m = (j-1)*(2*res+1) + i;
%  
%     
%    if pixel_intersections(m ,1) == 0
%        pixel_intersections(m, 2) = x;
%        pixel_intersections(m, 3) = y;
%    elseif pixel_intersections(m, 1) == 1
%        pixel_intersections(m, 4) = x;
%        pixel_intersections(m, 5) = y;       
%    end
%    pixel_intersections(m, 1) = pixel_intersections(m, 1) + 1;
%    
%    
%    if j < 2*res +1
%         if pixel_intersections(m+2*res+1+1 ,1) == 0
%             pixel_intersections(m+2*res+1+1, 2) = x;
%             pixel_intersections(m+2*res+1+1, 3) = y;
%         elseif pixel_intersections(m+2*res+1+1, 1) == 1
%             pixel_intersections(m+2*res+1+1, 4) = x;
%             pixel_intersections(m+2*res+1+1, 5) = y;       
%         end
%         pixel_intersections(m+2*res+1+1, 1) = pixel_intersections(m+2*res+1+1, 1) + 1;
%    end
% end
% 
% for n=1:M
%     if pixel_intersections(n, 1) == 2
%         x1 = pixel_intersections(n, 2);
%         y1 = pixel_intersections(n, 3);
%         x2 = pixel_intersections(n, 4);
%         y2 = pixel_intersections(n, 5);
% 
%         lengths(n) = sqrt((x1 - x2)^2+(y1 - y2)^2);
%     end
% end


end

