function O = and_(O,S,varargin)
% and_ - overloads '&' operator, computes the intersection of an empty set
%    and another set or numerical vector
%
% Syntax:
%    O = and_(O,S)
%    O = and_(O,S,method)
%
% Inputs:
%    O - emptySet object
%    S - contSet object or numerical vector
%    method - (optional) approximation method
%
% Outputs:
%    O - intersection
%
% Example: 
%    O = emptySet(2);
%    S = zonotope([1;1],[2 1; -3 1]);
%    O & S
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   05-April-2023 (MW, second argument can be empty)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% intersection is always empty set
[O,S] = findClassArg(O,S,'emptySet');

% ------------------------------ END OF CODE ------------------------------
