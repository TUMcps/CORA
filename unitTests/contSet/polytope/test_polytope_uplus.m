function res = test_polytope_uplus
% test_polytope_uplus - unit test function of uplus
%
% Syntax:
%    res = test_polytope_uplus
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

% Authors:       Tobias Ladner
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test empty case
assert(representsa(+polytope.empty(2),'emptySet'));

% 2D, bounded
A = [1 0 -1 0 1; 0 1 0 -1 1]'; b = [3; 2; 3; 2; 1];
P = polytope(A,b);
pP = +P;
% everything should remain as is
assert(all(pP.A == A, 'all'));
assert(all(pP.b == b, 'all'));
assert(isequal(pP, P));

% 2D, fully empty
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
pP = +P;
assert(representsa(pP,'fullspace'));


% add results
res = true;

% ------------------------------ END OF CODE ------------------------------
