function res = test_zonotope_underapproximate
% test_zonotope_underapproximate - unit test function of underapproximate
%
% Syntax:
%    res = test_zonotope_underapproximate
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
 
% assume true
res = true;

% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

% create direction matrix
S = [1 1; 0 1];

% obtain result 1
V_1 = underapproximate(Z1);

% obtain result 2
V_2 = underapproximate(Z1,S);

% true result 1
true_V_1 = [-10, 2, -10, 2; ...
            10, -8, 10, -8];
        
% true result 2
true_V_2 = [2, -10, -4, -4; ...
    -8, 10, 6, -4];

% check results
assert(compareMatrices(V_1,true_V_1) && compareMatrices(V_2,true_V_2));

% ------------------------------ END OF CODE ------------------------------
