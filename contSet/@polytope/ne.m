function res = ne(P,S,varargin)
% ne - Overloaded '~=' operator for the comparison of a polytope and
%    another contSet object
%
% Syntax:
%    res = P ~= S
%    res = ne(P,S)
%    res = ne(P,S,tol)
%
% Inputs:
%    P - polytope object 
%    S - polytope object 
%    tol - (optional) tolerance
%
% Outputs:
%    res - result of comparison
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([0 -1; 0 1;-1 0; 1 0; -1 -1],[2;3;2;3;2]);
%    res = P1 ~= P2;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       27-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call isequal (avoids duplicate code)
res = ~isequal(P,S,varargin{:});

% ------------------------------ END OF CODE ------------------------------
