function res = test_zonotope_rank
% test_zonotope_rank - unit test function of rank
%
% Syntax:
%    res = test_zonotope_rank
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

% Authors:       Matthias Althoff
% Written:       26-July-2016
% Last update:   15-September-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 2D zonotope
Z = zonotope([1, 2, 0, 4; 5, 6, 0, 0; -1, 4, 0, 8]);
res(end+1,1) = rank(Z) == 2;

% empty zonotope
Z_empty = zonotope.empty(2);
res(end+1,1) = rank(Z_empty) == 0;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
