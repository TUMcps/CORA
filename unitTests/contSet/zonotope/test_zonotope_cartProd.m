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

res = true(0);

% zonotope-zonotope case:

% 2D and 1D zonotopes
Z1 = zonotope([1,2,3,4; 5 6 7 8]);
Z2 = zonotope([9 10 11]);
Z_ = cartProd(Z1,Z2);
c = Z_.c; G = Z_.G;
% compare to true result
true_c = [1; 5; 9];
true_G = [2, 3, 4, 0, 0; ...
          6, 7, 8, 0, 0; ...
          0, 0, 0, 10,11];
res(end+1,1) = compareMatrices(c,true_c) && compareMatrices(G,true_G);


% 2D zonotope, 1D numeric
Z1 = zonotope([0;2],[3 4 2; -3 -1 3]);
num = 1;
Z_ = cartProd(Z1,num);
c = Z_.c; G = Z_.G;
% compare to true result
true_c = [0; 2; 1];
true_G = [3  4 2; ...
         -3 -1 3; ...
          0  0 0];
res(end+1,1) = compareMatrices(c,true_c) && compareMatrices(G,true_G);

% other ordering
Z_ = cartProd(num,Z1);
c = Z_.c; G = Z_.G;
% compare to true result
true_c = [1; 0; 2];
true_G = [0  0 0; ...
          3  4 2; ...
         -3 -1 3];
res(end+1,1) = compareMatrices(c,true_c) && compareMatrices(G,true_G);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
