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
% See also: -

% Authors:       Tobias Ladner
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = true(0);

% init
C = [1 0 -1 0 1; 0 1 0 -1 1]';
d = [3; 2; 3; 2; 1];
P = polytope(C,d);

% plus
pP = +P;
resvec(end+1) = all(pP.A == C, 'all');
resvec(end+1) = all(pP.b == d, 'all');

% compare with P
resvec(end+1) = isequal(pP, P);

% test empty case
resvec(end+1) = representsa(+polytope(),'emptySet');

% add results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
