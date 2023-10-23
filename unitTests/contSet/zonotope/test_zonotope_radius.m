function res = test_zonotope_radius
% test_zonotope_radius - unit test function of radius
%
% Syntax:
%    res = test_zonotope_radius
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
% Written:       27-August-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create a zonotope
Z = zonotope(zeros(2,1),[1 3; 2 1]);

% compute radius
r = radius(Z);

% analytical solution
r_true = 5;

% check results
res = withinTol(r,r_true);

% ------------------------------ END OF CODE ------------------------------
