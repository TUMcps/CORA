function res = test_emptySet_generateRandom
% test_emptySet_generateRandom - unit test function of generateRandom
%
% Syntax:
%    res = test_emptySet_generateRandom
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

% Authors:       Tobias Ladner
% Written:       02-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% without input arguments
O = emptySet.generateRandom();

% dimension given
n = 2;
O = emptySet.generateRandom('Dimension',n);
res = dim(O) == n;

% ------------------------------ END OF CODE ------------------------------
