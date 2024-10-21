function res = test_emptySet_isFullDim
% test_emptySet_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = test_emptySet_isFullDim
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
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% init empty set
n = 2;
O = emptySet(n);

% check result
assert(~isFullDim(O));

% ------------------------------ END OF CODE ------------------------------
