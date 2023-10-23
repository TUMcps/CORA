function res = test_zonotope_enclose
% test_zonotope_enclose - unit test function of enclose
%
% Syntax:
%    res = test_zonotope_enclose
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
% Last update:   09-August-2020 (MW, enhance randomness)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% create zonotopes
Z1 = zonotope([1,2,3,4; 5 6 7 8]);
Z2 = zonotope([9, 10, 11; 12, 13, 14]);

% obtain enclosing zonotope
Z_ = enclose(Z1,Z2);

% obtain zonotope matrix
c_ = Z_.c;
G_ = Z_.G;

% true result
true_c = [5; 8.5];
true_G = [6, 7, -4, -4, -4, 4; ...
            9.5, 10.5, -3.5, -3.5, -3.5, 8];

% check result
resvec(end+1) = compareMatrices(c_,true_c) && compareMatrices(G_,true_G);

% compute either
Z12 = enclose(Z1,Z2);
Z21 = enclose(Z1,Z2);

resvec(end+1) = isequal(Z12, Z21);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
