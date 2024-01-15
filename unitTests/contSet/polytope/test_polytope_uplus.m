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

resvec = true(0);

% test empty case
resvec(end+1) = representsa(+polytope.empty(2),'emptySet');

% 2D, bounded
A = [1 0 -1 0 1; 0 1 0 -1 1]'; b = [3; 2; 3; 2; 1];
P = polytope(A,b);
pP = +P;
% everything should remain as is
resvec(end+1) = all(pP.A == A, 'all');
resvec(end+1) = all(pP.b == b, 'all');
resvec(end+1) = isequal(pP, P);

% 2D, fully empty
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
pP = +P;
resvec(end+1) = representsa(pP,'fullspace');


% add results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
