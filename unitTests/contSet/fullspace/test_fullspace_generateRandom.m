function res = test_fullspace_generateRandom
% test_fullspace_generateRandom - unit test function of generateRandom
%
% Syntax:
%    res = test_fullspace_generateRandom
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
% Written:       25-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% without input arguments
fs = fullspace.generateRandom();

% dimension given
n = 2;
fs = fullspace.generateRandom('Dimension',n);
assert(dim(fs) == n);

% ------------------------------ END OF CODE ------------------------------
