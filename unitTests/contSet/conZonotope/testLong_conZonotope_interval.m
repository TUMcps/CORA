function res = testLong_conZonotope_interval
% testLong_conZonotope_interval - unit test function for the
%    calculation of a bounding box of a constrained zonotope object
%
% Syntax:
%    res = testLong_conZonotope_interval
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
% See also: -
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Niklas Kochdumper
% Written:       22-May-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Generate random conZonotope object
points = rand(3,100);
ind = convhulln(points');
ind = unique(ind(:,1),'stable');
V = points(:,ind);

P = polytope(V);
cZ = conZonotope(P);

% calculate interval
I = interval(cZ);

% compare with ground-truth for the vertices
V = vertices(cZ);
I_ = interval(min(V,[],2),max(V,[],2));

for i = 1:length(I)
    if ~isequal(I,I_,1e-10)
        res = false;
        break
    end
end

% ------------------------------ END OF CODE ------------------------------
