function cZ = conZonotope(zB)
% conZonotope - convert a zonotope bundle to a constrained zonotope
%
% Syntax:
%    cZ = conZonotope(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    cZ - conZonotope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       23-May-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if zB.parallelSets == 0
    cZ = conZonotope.empty(dim(zB)); return
end

% initialization
cZ = conZonotope(zB.Z{1});

% calculate the intersection of the parallel sets
for i = 2:zB.parallelSets
    cZ_i = conZonotope(zB.Z{i});
    cZ = and_(cZ,cZ_i,'exact');
end

% ------------------------------ END OF CODE ------------------------------
