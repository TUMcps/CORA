function [S1,S2] = reorder(S1,S2)
% reorder - reorder two sets according to their preference value; numeric
%    types are also reordered to position 2
%
% Syntax:
%    [S1,S2] = reorder(S1,S2)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object
%
% Outputs:
%    S1 - contSet object with lower precendence
%    S2 - contSet object with higher precendence
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       28-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% classic swap using temporary variable
if isnumeric(S1) || ...
        isa(S2,'contSet') && isscalar(S2) && S2.precedence < S1.precedence
    S1_copy = S1;
    S1 = S2;
    S2 = S1_copy;
end

% ------------------------------ END OF CODE ------------------------------
