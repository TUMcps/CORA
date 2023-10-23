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

% no input argument (default: dimension 1)
fs = fullspace();
res = fs.dimension == 0;

% dimension greater or equal to one
n = 2;
fs = fullspace(n);
res(end+1,1) = fs.dimension == n;

% too many input arguments
try
    fs = fullspace(n,n);
    res(end+1,1) = false;
catch
    res(end+1,1) = true;
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
