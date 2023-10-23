function res = test_fullspace_isemptyobject
% test_fullspace_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_fullspace_isemptyobject
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

% init fullspace
n = 2;
fs = fullspace(n);

% check emptiness
res = ~isemptyobject(fs);

% ------------------------------ END OF CODE ------------------------------
