function zB = project(zB,dims)
% project - projects a zonotope bundle onto the specified dimensions; note
%    that this returns an outer-approximation
%
% Syntax:
%    zB = project(zB,dims)
%
% Inputs:
%    zB - (zonoBundle) zonotope bundle
%    dims - dimensions for projection
%
% Outputs:
%    zB - (zonoBundle) projected zonotope bundle
%
% Example: 
%    Z1 = zonotope(zeros(3,1),[1 0.5 0; -0.2 1 0.3; -1 0.4 -0.2]);
%    Z2 = zonotope(ones(3,1),[1 -0.5 1; 0.2 1 -0.4; 0.5 -0.7 0.2]);
%    zB = zonoBundle({Z1,Z2});
%    res = project(zB,[1,2]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       04-February-2011
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% project for each zonotope
for i=1:zB.parallelSets
    zB.Z{i} = project(zB.Z{i},dims);
end

% ------------------------------ END OF CODE ------------------------------
