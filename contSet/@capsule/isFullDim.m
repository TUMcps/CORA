function res = isFullDim(C)
% isFullDim - checks if the dimension of the affine hull of a capsule is
%    equal to the dimension of its ambient space
%
% Syntax:
%    res = isFullDim(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    res - true/false
%
% Example:
%    C = capsule([2;1],[1;1],1);
%    isFullDim(C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       02-January-2020 
% Last update:   09-January-2024 (MW, simplify)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = ~isempty(C.r) && C.r > 0;

% ------------------------------ END OF CODE ------------------------------
