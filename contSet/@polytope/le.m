function res = le(P,S)
% le - overloads '<=' operator, checks if P is contained in S
%
% Syntax:
%    le(P1,P2)
%
% Inputs:
%    P - polytope object
%    S - contSet object
%
% Outputs:
%    P - polytope object
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[4;6;4;6;4]);
%
%    res = le(P1,P2)
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

res = contains_(S,P,'exact');

% ------------------------------ END OF CODE ------------------------------
