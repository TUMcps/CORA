function res = test_zonotope_supportFunc
% test_zonotope_supportFunc - unit test function of support function
%
% Syntax:
%    res = test_zonotope_supportFunc
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

% Authors:       Mark Wetzlinger, Victor Gassmann
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty set
Z_empty = zonotope.empty(2);
assert(supportFunc(Z_empty,[1;1],'upper') == -Inf);
assert(supportFunc(Z_empty,[1;1],'lower') == Inf);

% instantiate zonotope
c = [2; 1];
G = [4 2 -2; 1 3 7];
Z = zonotope(c,G);

% check a couple of evaluations
assert(withinTol(supportFunc(Z,[1;0]),10));
assert(withinTol(supportFunc(Z,[-1;0]),6));
assert(withinTol(supportFunc(Z,[0;1]),12));
assert(withinTol(supportFunc(Z,[0;-1]),10));

% check 'range'
assert(isequal(supportFunc(Z,[1;0],'range'),...
        interval(supportFunc(Z,[1;0],'lower'),supportFunc(Z,[1;0],'upper'))));

% check a couple of support vectors
[~,x1] = supportFunc(Z,[1;1]);
[~,x2] = supportFunc(Z,[-1;1]);
[~,x3] = supportFunc(Z,[-1;-1]);
assert(compareMatrices([x1 x2 x3],[6 12; -2 10; -2 -10]'));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
