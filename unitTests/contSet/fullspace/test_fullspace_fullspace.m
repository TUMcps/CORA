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

res = true(0);

% dimension greater or equal to one
n = 2;
fs = fullspace(n);
res(end+1,1) = fs.dimension == n;

% combine results
res = all(res);

% too many input arguments
try
    fs = fullspace(n,n);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
