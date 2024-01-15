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

res = true(0);

% empty capsule
C = capsule.empty(2);
res(end+1,1) = isempty(radius(C));

% 2D capsule
C = capsule([0;0],[1;0],0.5);
res(end+1,1) = withinTol(radius(C),1.5);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
