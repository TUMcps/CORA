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

res = true(0);

% empty interval
I = interval.empty(2);
res(end+1,1) = isemptyobject(I);

% 3D interval
lb = [-2;-1;-3]; ub = [1;1;2];
I = interval(lb,ub);
res(end+1,1) = ~isemptyobject(I);

% interval matrix
lb = [-2 0; -1 1; -3 -2]; ub = [1 1; 1 2; 2 0];
I = interval(lb,ub);
res(end+1,1) = ~isemptyobject(I);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
