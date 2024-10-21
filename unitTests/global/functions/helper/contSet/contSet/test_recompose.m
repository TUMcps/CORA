function res = test_recompose()
% test_recompose - unit test function for the recomposition of a cell array
%    of sets to a single set
%
% Syntax:
%    res = test_recompose()
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
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-12;

% intervals
S = {interval(1,2); ...
     interval([2;-1],[4;5]); ...
     interval(2,3); ...
     interval([-1;0;1],[2;0;2])};

S_rec = recompose(S);
S_true = interval([1;2;-1;2;-1;0;1], [2;4;5;3;2;0;2]);
assert(isequal(S_rec,S_true,tol));

% zonotopes
S = {zonotope(1,2); ...
     zonotope([2;-1],[4 0; -1 5]); ...
     zonotope(2,3)};

S_rec = recompose(S);
S_true = zonotope([1;2;-1;2], [2 0 0 0; 0 4 0 0; 0 -1 5 0; 0 0 0 3]);
assert(isequal(S_rec,S_true,tol));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
