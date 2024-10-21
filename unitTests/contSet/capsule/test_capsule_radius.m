function res = test_capsule_radius
% test_capsule_radius - unit test function of radius
%
% Syntax:
%    res = test_capsule_radius
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
% Written:       28-August-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty capsule
C = capsule.empty(2);
assert(isempty(radius(C)));

% 2D capsule
C = capsule([0;0],[1;0],0.5);
assert(withinTol(radius(C),1.5));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
