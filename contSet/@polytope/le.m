function res = le(P,S,varargin)
% le - overloads '<=' operator, checks if P is contained in S
%
% Syntax:
%    res = le(P1,P2)
%    res = le(P1,P2,tol)
%
% Inputs:
%    P - polytope object
%    S - contSet object
%    tol - (optional) tolerance
%
% Outputs:
%    P - polytope object
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[4;6;4;6;4]);
%
%    res = P1 <= P2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Viktor Kotsev
% Written:       09-May-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set tolerance
tol = setDefaultValues({eps},varargin);

% call contains and handle everything there
res = contains_(S,P,'exact',tol);

% ------------------------------ END OF CODE ------------------------------
