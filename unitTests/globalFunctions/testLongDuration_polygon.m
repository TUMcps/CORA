function res = testLongDuration_polygon()
% testLongDuration_polygon - unit test function for the polgon class
%
% Syntax:  
%    res = testLongDuration_polygon()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon

% Author:       Niklas Kochdumper
% Written:      16-June-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;


% Test 1: Quadratic Map ---------------------------------------------------

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
    if ~in(res,points_)
        
        % save variables so that failure can be reproduced
        path = pathFailedTests(mfilename());
        save(path,'pgon','Q','points');

        error('Random test failed!');
    end
end

res = true;

%------------- END OF CODE --------------