function res = test_fullspace_Inf
% test_fullspace_Inf - unit test function of R^n instantiation
%
% Syntax:
%    res = test_fullspace_Inf
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

% Authors:       Tobias Ladner
% Written:       16-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1D
n = 1;
fs = fullspace.Inf(n);
assert(representsa(fs,'fullspace') && dim(fs) == 1);

% 5D
n = 5;
fs = fullspace.Inf(n);
assert(representsa(fs,'fullspace') && dim(fs) == 5);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
