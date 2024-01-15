function res = test_capsule_display
% test_capsule_display - unit test function of display
%
% Syntax:
%    res = test_capsule_display
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
% Written:       24-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% empty capsule
n = 2;
C = capsule.empty(n)

% set center, generator, and radius
c = [2; 0];
g = [1; -1];
r = 0.5;

% only center
C = capsule(c)

% center and generator
C = capsule(c,g)

% center, generator, and radius
C = capsule(c,g,r)

% ------------------------------ END OF CODE ------------------------------
