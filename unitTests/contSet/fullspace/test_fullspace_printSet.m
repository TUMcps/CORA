function res = test_fullspace_printSet
% test_fullspace_printSet - unit test function of printSet
%
% Syntax:
%    res = test_fullspace_printSet
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

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test normal set
n = 2;
fs = fullspace(n);

printSet(fs)
printSet(fs,'high')
printSet(fs,'high',true)
printSet(fs,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
