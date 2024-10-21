function [A,b] = priv_equalityToInequality(A,b,Ae,be)
% priv_equalityToInequality - rewrites all equality constraints as
%    inequality constraints
%
% Syntax:
%    [A,b] = priv_equalityToInequality(A,b,Ae,be)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%
% Outputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

A = [A; Ae; -Ae];
b = [b; be; -be];

% ------------------------------ END OF CODE ------------------------------
