function S_out = plus(O,S)
% plus - overloaded '+' operator for the Minkowski addition of a
%    full-dimensional space and another set or vector
%
% Syntax:
%    S_out = O + S
%    S_out = plus(O,S)
%
% Inputs:
%    O - emptySet object, numeric
%    S - contSet object, numeric
%
% Outputs:
%    S_out - Minkowski sum
%
% Example: 
%    O = emptySet(2);
%    Z = zonotope([1;1],[2 1; -1 0]);
%    O + Z
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[S_out,S] = reorderNumeric(O,S);

% check dimensions of ambient space
equalDimCheck(S_out,S);

% ------------------------------ END OF CODE ------------------------------
