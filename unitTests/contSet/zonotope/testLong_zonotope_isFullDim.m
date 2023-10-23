function res = testLong_zonotope_isFullDim
% testLong_zonotope_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = testLong_zonotope_isFullDim
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
% Written:       12-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 1000;

for i=1:nrTests

    % random dimension
    n = randi([2,50]);

    % random center
    c = randn(n,1);

    % random number of generators
    gamma = randi([n,100]);

    % random full-rank generator matrix
    G = randn(n,gamma);
    while rank(G) < n
        G = randn(n,gamma);
    end
    
    % instantiate zonotope
    Z = zonotope(c,G);
    
    % check if full-dimensional
    if ~isFullDim(Z)
        res = false; return
    end
end

% ------------------------------ END OF CODE ------------------------------
