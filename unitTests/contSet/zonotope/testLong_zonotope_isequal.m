function res = testLong_zonotope_isequal
% testLong_zonotope_isequal - unit test function of isequal
%
% Syntax:
%    res = testLong_zonotope_isequal
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
% Written:       17-September-2019
% Last update:   21-April-2020
%                09-August-2020 (enhance randomness)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 1000;

% compare randomly generated zonotopes
for i=1:nrTests
    % random dimension
    n = randi(10);

    % create two random zonotopes
    nrOfGens = randi([2*n,5*n]);
    Z1 = zonotope(-1+2*rand(n,nrOfGens+1));
    Z2 = zonotope(-1+2*rand(n,nrOfGens+1));

    % check all combinations
    if ~(isequal(Z1,Z1) && isequal(Z2,Z2) && ~isequal(Z1,Z2))
        res = false; return
    end
    
    % check different order of generators
    c3 = zeros(n,1);
    G3 = -1+2*rand(n,nrOfGens);
    G3difforder = G3(:,randperm(nrOfGens));
    
    % init zonotopes
    Z3 = zonotope(c3,G3);
    Z3difforder = zonotope(c3,G3difforder);
    
    % compare zonotopes
    if ~isequal(Z3,Z3difforder)
        res = false; return
    end
end

% ------------------------------ END OF CODE ------------------------------
