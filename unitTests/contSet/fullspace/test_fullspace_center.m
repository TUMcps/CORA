function res = test_fullspace_center
% test_fullspace_center - unit test function of center
%
% Syntax:
%    res = test_fullspace_center
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
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init fullspace
n = 2;
fs = fullspace(n);

% compute center
c = center(fs);
c_true = zeros(n,1);

% compare results
res = all(c == c_true);

% ------------------------------ END OF CODE ------------------------------
