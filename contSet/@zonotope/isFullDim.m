function res = isFullDim(Z)
% isFullDim - checks if the dimension of the affine hull of a zonotope is
%    equal to the dimension of its ambient space
%
% Syntax:
%    res = isFullDim(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    res - true/false
%
% Example:
%    Z1 = zonotope([1 2 1;3 1 2]);
%    Z2 = zonotope([1 2 1;3 4 2]);
%
%    isFullDim(Z1)
%    isFullDim(Z2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: isempty

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       02-January-2020 
% Last update:   12-March-2021 (MW, add empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~representsa_(Z,'emptySet',eps)
    Zdim = dim(Z);
    Grank = rank(Z.G);
    res = Zdim == Grank;
else
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
