function res = test_polytope_display
% test_polytope_display - unit test function of display
%
% Syntax:
%    res = test_polytope_display
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       04-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty polytope
P = polytope.empty(2)

% only inequalities
A = [2 1; 1 2; -1 2; -2 1; 0 -2];
b = ones(5,1);
P = polytope(A,b)

% only equalities
Ae = [1 1; -1 1];
be = [2; 2];
P = polytope([],[],Ae,be)

% both inequalities and equalities
P = polytope(A,b,Ae,be)

% high-dimensional system
A = [eye(10);-eye(10)];
b = ones(20,1);
P = polytope(A,b)


% all displayed successfully
res = true;

% ------------------------------ END OF CODE ------------------------------
