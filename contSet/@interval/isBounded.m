function res = isBounded(I)
% isBounded - determines if an interval is bounded
%
% Syntax:
%    res = isBounded(I)
%
% Inputs:
%    I - interval object
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

% unbounded if any limit is -Inf or +Inf, otherwise bounded
if isemptyobject(I)
    res = true;
else
    res = any(I.inf == -Inf) || any(I.sup == Inf);
end

% ------------------------------ END OF CODE ------------------------------
