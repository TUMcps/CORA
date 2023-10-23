function res = testLong_zonotope_dim
% testLong_zonotope_dim - unit test function of dim
%
% Syntax:
%    res = testLong_zonotope_dim
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
% Written:       11-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% number of tests
nrOfTests = 1000;

for i=1:nrOfTests
    % random dimension
    n = randi([2,50]);
    % random center
    c = randn(n,1);
    % random generator matrix (also empty)
    if rand(1) < 0.05
        G = [];
    else
        G = 5*randn(n,randi(10));
    end
    
    % instantiate zonotope
    Z = zonotope(c,G);
    
    % get dimension
    Zdim = dim(Z);
    
    % assert correctness
    if Zdim ~= n
        res = false; break
    end
end

% ------------------------------ END OF CODE ------------------------------
