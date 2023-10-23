function res = test_taylm_generateRandom
% test_taylm_generateRandom - unit test function of generateRandom
%
% Syntax:
%    res = test_taylm_generateRandom
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
tay = taylm.generateRandom();

% dimension given
n = 2;
tay = taylm.generateRandom('Dimension',n);
res = dim(tay) == n;

% ------------------------------ END OF CODE ------------------------------
