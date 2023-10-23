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

% create zonotope
Z = zonotope([1, 2, 0, 4; 5, 6, 0, 0; -1, 4, 0, 8]);

% obtain zonotope without zeros
d = rank(Z);

% true result
true_val = 2;

% check result
res = withinTol(d,true_val) && rank(zonotope()) == 0;

% ------------------------------ END OF CODE ------------------------------
