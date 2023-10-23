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

res = true;

% empty input arguments
zB = zonoBundle.generateRandom();

% only dimension
n = 3;
zB = zonoBundle.generateRandom('Dimension',n);
res(end+1,1) = dim(zB) == n;

% only number of parallel sets
nrZonos = 5;
zB = zonoBundle.generateRandom('NrZonotopes',nrZonos);
res(end+1,1) = zB.parallelSets == nrZonos;

% dimension and number of parallel sets
zB = zonoBundle.generateRandom('Dimension',n,'NrZonotopes',nrZonos);
res(end+1,1) = dim(zB) == n && zB.parallelSets == nrZonos;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
