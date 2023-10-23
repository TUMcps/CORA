function res = test_interval_isemptyobject
% test_interval_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_interval_isemptyobject
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

% Authors:       Mark Wetzlinger
% Written:       03-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate intervals
I1 = interval();

lb = [-2;-1;-3];
ub = [1;1;2];
I2 = interval(lb,ub);

% check results
res = isemptyobject(I1) && ~isemptyobject(I2);

% ------------------------------ END OF CODE ------------------------------
