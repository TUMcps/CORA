function res = test_capsule_randPoint
% test_capsule_randPoint - unit test function of randPoint
%
% Syntax:
%    res = test_capsule_randPoint
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required:
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty case
C = capsule.empty(2);
p = randPoint(C);
res(end+1,1) = isnumeric(p) && isempty(p) && all(size(p) == [2, 0]);

% degenerate capsule
c = [1; -1; 2]; g = [0; 0; 0]; r = 0;
C = capsule(c,g,r);
% rand point can only be the center
p = randPoint(C);
res(end+1,1) = compareMatrices(p,c);

% instantiate full-dimensional capsule
c = [3; 0; 0]; g = [-1; 3; 2]; r = 1;
C = capsule(c,g,r);
% generate random points
p = randPoint(C,10);
% check if all random points inside capsule
res(end+1,1) = all(contains(C,p));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
