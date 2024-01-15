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

% test empty case
resvec(end+1) = representsa(-polytope.empty(2),'emptySet');

% 2D, bounded
A = [1 0 -1 0 1; 0 1 0 -1 1]'; b = [3; 2; 3; 2; 1];
P = polytope(A,b);
% negate
nP = -P;
resvec(end+1) = all(nP.A == -A, 'all');
resvec(end+1) = all(nP.b == b, 'all');
% compare with -1 * P
resvec(end+1) = isequal(nP, -1*P);

% 2D, fully empty
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
nP = -P;
resvec(end+1) = representsa(nP,'fullspace');


% add results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
