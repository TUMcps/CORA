function res = test_polytope_uminus
% test_polytope_uminus - unit test function of uminus
%
% Syntax:
%    res = test_polytope_uminus
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
% See also: -

% Authors:       Tobias Ladner
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test empty case
assert(representsa(-polytope.empty(2),'emptySet'));

% 2D, bounded
A = [1 0 -1 0 1; 0 1 0 -1 1]'; b = [3; 2; 3; 2; 1];
P = polytope(A,b);
% negate
nP = -P;
assert(all(nP.A == -A, 'all'));
assert(all(nP.b == b, 'all'));
% compare with -1 * P
assert(isequal(nP, -1*P));

% 2D, fully empty
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
nP = -P;
assert(representsa(nP,'fullspace'));


% add results
res = true;

% ------------------------------ END OF CODE ------------------------------
