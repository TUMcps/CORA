function res = test_zonoBundle_generateRandom
% test_zonoBundle_generateRandom - unit test function of generateRandom
%
% Syntax:
%    res = test_zonoBundle_generateRandom
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
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty input arguments
zB = zonoBundle.generateRandom();

% 1D
n = 1;
zB = zonoBundle.generateRandom('Dimension',n);
assert(dim(zB) == n);

% only dimension
n = 3;
zB = zonoBundle.generateRandom('Dimension',n);
assert(dim(zB) == n);

% only number of parallel sets
nrZonos = 5;
zB = zonoBundle.generateRandom('NrZonotopes',nrZonos);
assert(zB.parallelSets == nrZonos);

% dimension and number of parallel sets
n = 3; nrZonos = 5;
zB = zonoBundle.generateRandom('Dimension',n,'NrZonotopes',nrZonos);
assert(dim(zB) == n && zB.parallelSets == nrZonos);


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
