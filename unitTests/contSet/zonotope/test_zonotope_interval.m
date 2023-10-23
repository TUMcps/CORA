function res = test_zonotope_interval
% test_zonotope_interval - unit test function of interval
%
% Syntax:
%    res = test_zonotope_interval
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

% create parallelotopes
I1 = interval(Z1);

% obtain results
lb = infimum(I1);
ub = supremum(I1);

% true results
true_lb = [-10; -8];
true_ub = [2; 10];   

% check result
res = compareMatrices(lb,true_lb) && compareMatrices(ub,true_ub);

% ------------------------------ END OF CODE ------------------------------
