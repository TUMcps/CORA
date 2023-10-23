function res = ne(P1,P2,varargin)
% ne - Overloaded '~=' operator for the comparison of polytopes
%
% Syntax:
%    res = P1 ~= P2
%    res = ne(P1,P2)
%    res = ne(P1,P2,tol)
%
% Inputs:
%    P1 - polytope object 
%    P2 - polytope object 
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
res = ~isequal(P1,P2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
