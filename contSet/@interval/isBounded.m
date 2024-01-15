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

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       24-July-2023
% Last update:   26-October-2023 (TL, reworked non-empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% unbounded if any limit is -Inf or +Inf, otherwise bounded
if representsa_(I,'emptySet',0)
    res = true;
else
    res = ~any(isinf(I.inf)) && ~any(isinf(I.sup));
end

% ------------------------------ END OF CODE ------------------------------
