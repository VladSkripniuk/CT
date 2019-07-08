function distances = CalculateDistances(pointlist, M)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

distances = zeros(M, 1);
res = round((sqrt(M) - 1) /2);

for k=1:size(pointlist, 1)-1
    x = (pointlist(k, 1)+pointlist(k+1, 1))/2;
    y = (pointlist(k, 2)+pointlist(k+1, 2))/2;
   [i, j] = floats2pixels(x, y, res);
   
   m = (j-1) * sqrt(M) + i;
   
%    if sqrt((pointlist(k,1)-pointlist(k+1, 1))^2 + (pointlist(k,2)-pointlist(k+1, 2))^2) < 1e-2
%        continue
%    end
%    
   distances(m) = sqrt((pointlist(k,1)-pointlist(k+1, 1))^2 + (pointlist(k,2)-pointlist(k+1, 2))^2);
end

end

