function res = test_conZonotope_supportFunc
% test_conZonotope_supportFunc - unit test function for the evaluation of
%    the support function
%
% Syntax:
%    res = test_conZonotope_supportFunc
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
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Mark Wetzlinger
% Written:       10-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty set
cZ = conZonotope.empty(2);
assert(supportFunc(cZ,[1;1]) == -Inf);
assert(supportFunc(cZ,[1;1],'lower') == Inf);

% TEST 1: Figure 1 in [1] -------------------------------------------------

% construct constrained zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1];
b = 1;
cZ = conZonotope(Z,A,b);

% check a couple of evaluations
assert(withinTol(supportFunc(cZ,[1;0]),3.5) ...
        && withinTol(supportFunc(cZ,[0;1]),2.5) ...
        && withinTol(supportFunc(cZ,[1;0],'lower'),-2.5,1e-6) ...
        && withinTol(supportFunc(cZ,[0;1],'lower'),-1.5),1e-6);
% check range
assert(isequal(supportFunc(cZ,[1;0],'range'),interval(-2.5,3.5),1e-6));

% check a couple of support vectors
[~,x1] = supportFunc(cZ,[-0.2;1]);
[~,x2] = supportFunc(cZ,[-1;-1]);
[~,x3] = supportFunc(cZ,[1;-0.2]);
assert(compareMatrices([x1,x2,x3],[-0.5 2.5; -2.5 -1.5; 3.5 -0.5]'),1e-8);

% check direction with special handling
dir = [1;1];
a = dir' * cZ.G;
A = [A; a]; b = [b; 0];
cZ = conZonotope(Z,A,b);

% evaluate support function
val = supportFunc(cZ,dir);
assert(withinTol(val,b(2)));


% TEST 2: Figure 2 in [1] -------------------------------------------------

% construct constrained zonotope
Z = [0 1 0 1;0 1 2 -1];
A = [-2 1 -1];
b = 2;
cZ = conZonotope(Z,A,b);

% check a couple of evaluations
assert(withinTol(supportFunc(cZ,[1;0]),0) ...
        && withinTol(supportFunc(cZ,[0;1]),3,1e-8) ...
        && withinTol(supportFunc(cZ,[1;0],'lower'),-2) ...
        && withinTol(supportFunc(cZ,[0;1],'lower'),-2));

% check a couple of support vectors
[~,x1] = supportFunc(cZ,[-1;-1]);
[~,x2] = supportFunc(cZ,[1;0]);
[~,x3] = supportFunc(cZ,[-0;1]);
assert(compareMatrices([x1,x2,x3],[-2 -2; 0 0; -1 3]',1e-8));


% TEST 3 ------------------------------------------------------------------

% construct constrained zonotope
Z = [0 3 0 1 -2;0 0 2 1 1];
A = [0 0 0 1];
b = 0.5;
cZ = conZonotope(Z,A,b);

% check a couple of evaluations
assert(withinTol(supportFunc(cZ,[1;0]),3) ...
        && withinTol(supportFunc(cZ,[0;1]),3.5) ...
        && withinTol(supportFunc(cZ,[1;0],'lower'),-5) ...
        && withinTol(supportFunc(cZ,[0;1],'lower'),-2.5));

% check a couple of support vectors
[~,x1] = supportFunc(cZ,[1;1]);
[~,x2] = supportFunc(cZ,[0.2;-1]);
[~,x3] = supportFunc(cZ,[-0.2;1]);
assert(compareMatrices([x1,x2,x3],[3 3.5; 1 -2.5; -3 3.5]'));


% TEST 4 ------------------------------------------------------------------

% degenerate constrained zonotope
c = ones(2,1);
G = [1 0.5 -2; 0 0 0];
A = [-2 1 -1];
b = 2;
cZ = conZonotope(c,G,A,b);

% direction normal to a direction with zero extension
dir = [0;1];

% compute support function
[~,x] = supportFunc(cZ,dir);

% x2 direction not changed by generators
assert(withinTol(x(2),c(2)));


% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
