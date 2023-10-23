function res = isFullDim(zB)
% isFullDim - checks if the dimension of the affine hull of a zonotope
%    bundle is equal to the dimension of its ambient space
%
% Syntax:
%    res = isFullDim(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    res - true/false
%
% Example:
%    I1 = interval([0;0],[2;2]);
%    I2 = interval([0;1],[2;3]);
%    I3 = interval([0;2],[2;4]);
%
%    zB1 = zonoBundle({zonotope(I1),zonotope(I2)});
%    zB2 = zonoBundle({zonotope(I1),zonotope(I3)});
%
%    isFullDim(zB1)
%    isFullDim(zB2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Authors:       Niklas Kochdumper
% Written:       02-January-2020 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if zB.parallelSets == 0
    res = false; return
end

% loop over all zonotopes
for i = 1:zB.parallelSets
   if ~isFullDim(zB.Z{i})
      res = false; return
   end
end

% convert to polytope and call polytope/isFullDim
res = isFullDim(polytope(zB));

% ------------------------------ END OF CODE ------------------------------
