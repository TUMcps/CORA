function res = test_zonotope_cartProd
% test_zonotope_cartProd - unit test function of cartesian product
%
% Syntax:
%    res = test_zonotope_cartProd
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-July-2016
% Last update:   03-January-2023 (MW, add zonotope-numeric cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% zonotope-zonotope case:

% 2D and 1D zonotopes
Z1 = zonotope([1,2,3,4; 5 6 7 8]);
Z2 = zonotope([9 10 11]);
Z_cartProd = cartProd(Z1,Z2);
% compare to true result
c_true = [1; 5; 9];
G_true = [2, 3, 4, 0, 0; ...
          6, 7, 8, 0, 0; ...
          0, 0, 0, 10,11];
assert(compareMatrices(Z_cartProd.c,c_true));
assert(compareMatrices(Z_cartProd.G,G_true));


% 2D zonotope, 1D numeric
Z1 = zonotope([0;2],[3 4 2; -3 -1 3]);
num = 1;
Z_cartProd = cartProd(Z1,num);
% compare to true result
c_true = [0; 2; 1];
G_true = [3  4 2; ...
         -3 -1 3; ...
          0  0 0];
assert(compareMatrices(Z_cartProd.c,c_true));
assert(compareMatrices(Z_cartProd.G,G_true));

% other ordering
Z_cartProd = cartProd(num,Z1);
% compare to true result
c_true = [1; 0; 2];
G_true = [0  0 0; ...
          3  4 2; ...
         -3 -1 3];
assert(compareMatrices(Z_cartProd.c,c_true));
assert(compareMatrices(Z_cartProd.G,G_true));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
