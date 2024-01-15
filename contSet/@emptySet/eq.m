function res = eq(O,S,varargin)
% eq - overloads '=='-operator for empty sets
%
% Syntax:
%    res = O == S
%    res = eq(O,S)
%    res = eq(O,S,tol)
%
% Inputs:
%    O - emptySet object
%    S - contSet object or numerical vector
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    O = emptySet(2);
%    S = zonotope.empty(2);
%    res1 = O == O;
%    res2 = O == S;
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

res = isequal(O,S,varargin);

% ------------------------------ END OF CODE ------------------------------
