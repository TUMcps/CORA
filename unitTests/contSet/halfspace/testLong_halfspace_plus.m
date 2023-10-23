function res = testLong_halfspace_plus
% testLong_halfspace_plus - unit test function of plus
%
% Syntax:
%    res = testLong_halfspace_plus
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
% Written:       16-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Random tests
res = true;
nrTests = 1000;

for i=1:nrTests
    % random dimension
    n = randi(50);
    
    % random normal vector (unit length)
    c = randn(n,1);
    c = c / vecnorm(c,2);
    
    % random distance
    d = randn(1);
    
    % init halfspace
    h = halfspace(c,d);
    
    % random vector
    v = randn(n,1);
    
    % compute result
    h_plus = h + v;
    
    % true result
    d_true = d + c.' * v;
    h_true = halfspace(c,d_true);
    
    % compare results
    if ~isequal(h_plus,h_true)
        res = false; break;
    end       
    
end

% ------------------------------ END OF CODE ------------------------------
