function res = isBounded(hyp)
% isBounded - determines if a constrained hyperplane is bounded
%
% Syntax:
%    res = isBounded(hyp)
%
% Inputs:
%    hyp - conHyperplane object
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       24-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1D always unbounded
if dim(hyp) == 1
    res = true;
else
    % convert to polytope and check
    res = isBounded(polytope(hyp));
end

% ------------------------------ END OF CODE ------------------------------
