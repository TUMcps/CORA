function res = isFullDim(E)
% isFullDim - checks if the dimension of the affine hull of an ellipsoid is
%    equal to the dimension of its ambient space
%
% Syntax:
%    res = isFullDim(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    res - true/false
%
% Example:
%    E1 = ellipsoid(eye(2));
%    isFullDim(E1)
%
%    E2 = ellipsoid([1 0; 0 0]);
%    isFullDim(E2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       02-January-2020 
% Last update:   24-March-2022 (remove dependency on object property)
%                06-July-2022 (VG, support class array case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

inputArgsCheck({{E,'att','ellipsoid'}}); 

if representsa_(E,'emptySet',eps)
    res = false;
else
    res = rank(E) == dim(E);
end

% ------------------------------ END OF CODE ------------------------------
