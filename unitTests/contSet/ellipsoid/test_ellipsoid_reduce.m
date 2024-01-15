function res = test_ellipsoid_reduce
% test_ellipsoid_reduce - unit test function of reduce
%
% Syntax:
%    res = test_ellipsoid_reduce
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
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
E = ellipsoid.empty(2);
% reduce
E_ = reduce(E);
% should remain the same...
res = isequal(E,E_);

% init column interval
E = ellipsoid([2 0; 0 1],[1;-1]);
% reduce
E_ = reduce(E);
% check result
res(end+1,1) = isequal(E,E_);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
