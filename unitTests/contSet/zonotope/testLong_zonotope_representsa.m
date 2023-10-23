function res = testLong_zonotope_representsa
% testLong_zonotope_representsa - unit test function of representsa
%
% Syntax:
%    res = testLong_zonotope_representsa
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
% Written:       09-August-2020
% Last update:   ---
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
    
    c = zeros(n,1);
    GnonInt = -1+2*rand(n,nrOfGens);
    ZnonInt = zonotope(c,GnonInt);
    
    % generate random matrix with only one non-negative entry per column
    GInt = zeros(n,nrOfGens);
    % non-negative row for each generator
    idx = randi([1,n],1,nrOfGens);
    % linear indexing
    idx = idx + [0:n:n*(nrOfGens-1)];
    % random values
    vals = -1+2*rand(1,nrOfGens);
    % write random values in matrix
    GInt(idx) = vals;
    % init zonotope
    ZInt = zonotope(c,GInt);

    % check zonotopes
    if ~(representsa(ZInt,'interval') && (n == 1 || ~representsa(ZnonInt,'interval')))
        res = false; return
    end
    
end

% ------------------------------ END OF CODE ------------------------------
