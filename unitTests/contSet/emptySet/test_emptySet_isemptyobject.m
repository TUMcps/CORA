function res = test_emptySet_isemptyobject
% test_emptySet_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_emptySet_isemptyobject
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
% Written:       25-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init empty set
n = 2;
O = emptySet(n);

% check result
res = ~isemptyobject(O);

% ------------------------------ END OF CODE ------------------------------
