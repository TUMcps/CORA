function res = test_conZonotope_vertices
% test_conZonotope_vertices - unit test function for the caclulation of
%                             vertices of a constrained zonotope object
%
% Syntax:
%    res = test_conZonotope_vertices
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
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Niklas Kochdumper
% Written:       11-May-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty set
n = 2;
cZ = conZonotope.empty(n);
V = vertices(cZ);
res(end+1,1) = isnumeric(V) && isempty(V) && size(V,1) == n;

% TEST 1: Figure 1 in [1] -------------------------------------------------

% construct zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1];
b = 1;
cZ = conZonotope(Z,A,b);

% calculate vertices
V = vertices(cZ);

% plot the result
% plotZono(cZono);

% compare with ground-truth
V_ = [-0.5 3.5 -2.5;2.5 -0.5 -1.5];
res(end+1,1) = compareMatrices(V,V_);


% TEST 2: Figure 2 in [1] -------------------------------------------------

% construct zonotope
Z = [0 1 0 1;0 1 2 -1];
A = [-2 1 -1];
b = 2;
cZ = conZonotope(Z,A,b);

% calculate vertices
V = vertices(cZ);

% plot the result
% plotZono(cZono);

% compare with ground-truth
V_ = [-1 0 -2;3 0 -2];
res(end+1,1) = compareMatrices(V,V_);


% TEST 3 ------------------------------------------------------------------

% construct zonotope
Z = [0 3 0 1 -2;0 0 2 1 1];
A = [0 0 0 1];
b = 0.5;
cZ = conZonotope(Z,A,b);

% calculate vertices
V = vertices(cZ);

% plot the result
% plotZono(cZono);

% compare with ground-truth
V_ = [-4 -4 -2 4 4 2;-3 1 3 3 -1 -3] + [-1;0.5]*ones(1,6);
res(end+1,1) = compareMatrices(V,V_);


% TEST 4 ------------------------------------------------------------------

% construct zonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZ = conZonotope(Z,A,b);

% calculate vertices
V = vertices(cZ);

% plot the result
% plotZono(cZono);

% compare with ground-truth
V_ = [1 1 3 3;3 -1 2 -2];
res(end+1,1) = compareMatrices(V,V_);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
