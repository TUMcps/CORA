function res = test_zonotope_isequal
% test_zonotope_isequal - unit test function of isequal
%
% Syntax:
%    res = test_zonotope_isequal
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
% Written:       17-September-2019
% Last update:   21-April-2020
%                09-August-2020 (enhance randomness)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% create zonotopes
Z1 = zonotope([1, 2, 4;
               5, 6, 0;
              -1, 4, 8]);
Z2 = zonotope([1, 2, 4;
               4, 6, 0;
              -1, 4, 8]);
Z3 = zonotope([1, 2, 0, 4;
               5, 6, 0, 0;
              -1, 4, 0, 8]);

% check result
res(end+1,1) = isequal(Z1,Z3);
res(end+1,1) = ~isequal(Z1,Z2);

% different order of generators
Z1 = zonotope(ones(2,1),[1 2 5 3 3;
                          2 3 0 4 1]);
Z2 = zonotope(ones(2,1),[2 1 3 5 3;
                          3 2 4 0 1]);

res(end+1,1) = isequal(Z1,Z2);

% different sign
Z1 = zonotope(0,1);
Z2 = zonotope(0,-1);
res(end+1,1) = isequal(Z1,Z2);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
