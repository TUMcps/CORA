function res = testLong_polygon_quadMap()
% testLong_polygon_quadMap - unit test function for polygon/quadMap
%
% Syntax:
%    res = testLong_polygon_quadMap()
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

% Authors:       Niklas Kochdumper
% Written:       16-June-2021
% Last update:   ---
% Last revision: 25-May-2023 (TL, split unit tests)

% ------------------------------ BEGIN CODE -------------------------------

% test if the function quadMap is over-approximative

for i = 1:3
   
    % create random polygon and matrices for quadratic map
    pgon = polygon.generateRandom();
    
    Q = cell(2,1);
    Q{1} = -1 + 2*rand(2);
    Q{2} = -1 + 2*rand(2);
    
    % compute quadratic map
    res = quadMap(pgon,Q);
    
    % compute quadratic map for random points inside the polygon
    points = [randPoint(pgon,'all','extreme'),randPoint(pgon,1000)];
    points_ = zeros(size(points));
    
    for j = 1:size(points,2)
        points_(1,j) = points(:,j)'*Q{1}*points(:,j);
        points_(2,j) = points(:,j)'*Q{2}*points(:,j);
    end
    
    % check correctness
    if ~contains(res,points_)
        throw(CORAerror('CORA:testFailed'));
    end
end

res = true;

% ------------------------------ END OF CODE ------------------------------
