function res = test_fullspace_fullspace
% test_fullspace_fullspace - unit test function of constructor
%
% Syntax:
%    res = test_fullspace_fullspace
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

% dimension greater or equal to one
n = 2;
fs = fullspace(n);
assert(fs.dimension == n);

% too many input arguments
assertThrowsAs(@fullspace,'CORA:numInputArgsConstructor',n,n);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
