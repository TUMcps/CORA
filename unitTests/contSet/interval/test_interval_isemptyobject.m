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

% empty interval
I = interval.empty(2);
assert(isemptyobject(I));

% 3D interval
lb = [-2;-1;-3]; ub = [1;1;2];
I = interval(lb,ub);
assert(~isemptyobject(I));

% interval matrix
lb = [-2 0; -1 1; -3 -2]; ub = [1 1; 1 2; 2 0];
I = interval(lb,ub);
assert(~isemptyobject(I));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
