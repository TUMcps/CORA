function res = plus(summand1,summand2)
% plus - overloaded '+' operator for the Minkowski addition of a
%    full-dimensional space and another set or vector
%
% Syntax:
%    fs = plus(summand1,summand2)
%
% Inputs:
%    summand1 - emptySet object, contSet object, numerical vector
%    summand2 - emptySet object, contSet object, numerical vector
%
% Outputs:
%    res - emptySet object
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

% determine emptySet object
[res,summand] = findClassArg(summand1,summand2,'emptySet');

% check dimensions of ambient space
equalDimCheck(res,summand);

% ------------------------------ END OF CODE ------------------------------
