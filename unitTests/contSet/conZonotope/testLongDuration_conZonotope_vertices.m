function res = testLongDuration_conZonotope_vertices
% testLongDuration_conZonotope_vertices - unit test function for the 
%    calculation of vertices of a constrained zonotope object
%
% Syntax:  
%    res = testLongDuration_conZonotope_vertices
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Niklas Kochdumper
% Written:      24-August-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% loop over all test cases
for i = 1:100
    
    % generate random 2-dimensional constrained zonotope
    cZ = conZonotope.generateRandom('Dimension',2);

    % compute vertices
    V = vertices(cZ);

    % compute random extreme points of the constrained zonotope
    p = randPoint(cZ,100,'extreme');

    % check if all points are in the convex hull of the vertices
    if ~contains(mptPolytope(V'), p)
        path = pathFailedTests(mfilename());
        save(path,'cZ');
        res = false;
    end
end

%------------- END OF CODE --------------