function res = isnan(hyp)
% isnan - checks if a constrained hyperplane is defined using NaN
%
% Syntax:
%    res = isnan(hyp)
%
% Inputs:
%    hyp - conHyperplane object
%
% Outputs:
%    res - false
%
% Example: 
%    hyp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    res = isnan(hyp)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       01-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% NaN values are not possible by constructor
res = false;

% ------------------------------ END OF CODE ------------------------------
