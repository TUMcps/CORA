function res = isFullDim(I)
% isFullDim - checks if the dimension of the affine hull of an interval is
%    equal to the dimension of its ambient space
%
% Syntax:
%    res = isFullDim(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - true/false
%
% Example:
%    I1 = interval([-1;-2],[1;2]);
%    I2 = interval([-1;-2],[1;-2]);
%
%    isFullDim(I1)
%    isFullDim(I2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       02-January-2020 
% Last update:   12-March-2021 (MW, empty interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if representsa_(I,'emptySet',eps)
    res = false;
else
    res = ~any(withinTol(rad(I),0)); 
end

% ------------------------------ END OF CODE ------------------------------
