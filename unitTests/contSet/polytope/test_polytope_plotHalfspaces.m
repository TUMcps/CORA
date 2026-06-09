function res = test_polytope_plotHalfspaces
% test_polytope_plotHalfspaces - unit test function of plotHalfspaces
%
% Syntax:
%    res = test_polytope_plotHalfspaces
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
% Written:       22-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init polytope
A = [1 0 -1 0 1; 0 1 0 -1 1]';
b = [3; 2; 3; 2; 1];
P = polytope(A,b);

figure; hold on;
ax = gca();
plot(P); enlargeAxis;
plotHalfspaces(P);

assert(numel(allchild(ax)) == size(A,1) + 1)

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
