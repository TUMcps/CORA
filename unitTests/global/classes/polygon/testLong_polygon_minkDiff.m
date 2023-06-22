function res = testLong_polygon_minkDiff()
% testLong_polygon_minkDiff - unit test function for polygon/minkDiff
%
% Syntax:  
%    res = testLong_polygon_minkDiff()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon

% Author:       Niklas Kochdumper
% Written:      16-June-2021
% Last update:  ---
% Last revision:25-May-2023 (TL, split unit tests)

%------------- BEGIN CODE --------------

% test if the Minkowski difference yields a correct result

for i = 1:3
   
    % create random polygons
    pgon1 = polygon.generateRandom();
    pgon2 = polygon.generateRandom();
    
    scale = 0.1*rad(interval(pgon1))./rad(interval(pgon2));
    c = center(pgon2);
    pgon2 = c + diag(scale)*(pgon2 - c);
    
    % compute Minkowski difference
    pgon = minkDiff(pgon1,pgon2);
    
    % check correctness for random points
    points = randPoint(pgon1,100);
    
    for j = 1:size(points,2)
       if contains(pgon,points(:,j) - c)
          if ~contains(pgon1,points(:,j) - c + pgon2)
             throw(CORAerror('CORA:testFailed'));
          end
       else
          if contains(pgon2,points(:,j) - c + pgon2)
             throw(CORAerror('CORA:testFailed'));
          end
       end
    end 
end

res = true;

%------------- END OF CODE --------------