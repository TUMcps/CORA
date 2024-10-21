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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       17-September-2019
% Last update:   21-April-2020
%                09-August-2020 (enhance randomness)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% number of tests
nrTests = 1000;
tol = 1e-14;

% compare randomly generated zonotopes
for i=1:nrTests
    % random dimension
    n = randi(10);

    % create two random zonotopes
    nrOfGens = randi([2*n,5*n]);
    Z1 = zonotope(-1+2*rand(n,nrOfGens+1));
    Z2 = zonotope(-1+2*rand(n,nrOfGens+1));

    % check all combinations
    assertLoop(isequal(Z1,Z1,tol),i)
    assertLoop(isequal(Z2,Z2,tol),i)
    assertLoop(~isequal(Z1,Z2,tol),i)
    
    % check different order of generators
    c3 = zeros(n,1);
    G3 = -1+2*rand(n,nrOfGens);
    G3difforder = G3(:,randperm(nrOfGens));
    
    % init zonotopes
    Z3 = zonotope(c3,G3);
    Z3difforder = zonotope(c3,G3difforder);
    
    % compare zonotopes
    assertLoop(isequal(Z3,Z3difforder,tol),i)
end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
