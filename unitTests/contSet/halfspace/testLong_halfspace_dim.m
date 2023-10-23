function res = testLong_halfspace_dim
% testLong_halfspace_dim - unit test function of dim
%
% Syntax:
%    res = testLong_halfspace_dim
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
% Written:       27-September-2019
% Last update:   16-March-2021 (MW, add empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

nrTests = 1000;
res = true;

for i=1:nrTests
    % random dimension
    n = randi(50);
    
    % random normal vector and distance
    c = randn(n,1);
    d = randn(1);
    
    % init halfspace
    h = halfspace(c,d);
    
    % check result
    if dim(h) ~= n
        res = false; break;
    end
end

% ------------------------------ END OF CODE ------------------------------
