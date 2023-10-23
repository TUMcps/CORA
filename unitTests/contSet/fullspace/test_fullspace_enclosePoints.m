function res = test_fullspace_enclosePoints
% test_fullspace_enclosePoints - unit test function of enclosePoints
%
% Syntax:
%    res = test_fullspace_enclosePoints
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

% points
p = [2 -4 1; -4 3 2; 0 1 9];

% enclose by R^n
fs = fullspace.enclosePoints(p);

% check result
res = dim(fs) == 3;

% ------------------------------ END OF CODE ------------------------------
