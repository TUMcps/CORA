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

resvec = true(0);

% init
C = [1 0 -1 0 1; 0 1 0 -1 1]';
d = [3; 2; 3; 2; 1];
P = polytope(C,d);

% negate
nP = -P;
resvec(end+1) = all(nP.A == -C, 'all');
resvec(end+1) = all(nP.b == d, 'all');

% compare with -1 * P
resvec(end+1) = isequal(nP, -1*P);

% test empty case
resvec(end+1) = representsa(-polytope(),'emptySet');

% add results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
