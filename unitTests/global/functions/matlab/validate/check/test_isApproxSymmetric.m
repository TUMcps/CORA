function res = test_isApproxSymmetric
% test_isApproxSymmetric - unit test function for symmetry check
%
% Syntax:
%    res = test_isApproxSymmetric()
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
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% identity matrix
Q = eye(2);
assert(isApproxSymmetric(Q,0));

% slightly skewed
skew = 1e-10;
Q = [1 0; 0+skew, 0];
assert(~isApproxSymmetric(Q,0));
assert(isApproxSymmetric(Q,skew));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
