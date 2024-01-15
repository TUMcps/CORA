function res = eq(P,S,varargin)
% eq - Overloaded '==' operator for the comparison of polytopes
%
% Syntax:
%    res = P == S
%    res = eq(P,S)
%    res = eq(P,S,tol)
%
% Inputs:
%    P1 - polytope object 
%    P2 - contSet object or numerical vector
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([0 -1; 0 1;-1 0; 1 0; -1 -1],[2;3;2;3;2]);
%    res = P1==P2;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Viktor Kotsev
% Written:       13-June-2022 
% Last update:   27-April-2023 (MW, move code to isequal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if P1 and P2 contain each other
res = isequal(P,S,varargin{:});

% ------------------------------ END OF CODE ------------------------------
