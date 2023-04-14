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

% Author:       Niklas Kochdumper
% Written:      11-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% empty set
vertEmptycZ = vertices(conZonotope());
if ~isnumeric(vertEmptycZ) || ~isempty(vertEmptycZ)
    res = false;
end

% TEST 1: Figure 1 in [1] -------------------------------------------------

% construct zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1];
b = 1;
cZono = conZonotope(Z,A,b);

% calculate vertices
V = vertices(cZono);

% plot the result
% plotZono(cZono);

% compare with ground-truth
V_ = [-0.5 3.5 -2.5;2.5 -0.5 -1.5];

if ~compareMatrices(V,V_)
    res = false;
end



% TEST 2: Figure 2 in [1] -------------------------------------------------

% construct zonotope
Z = [0 1 0 1;0 1 2 -1];
A = [-2 1 -1];
b = 2;
cZono = conZonotope(Z,A,b);

% calculate vertices
V = vertices(cZono);

% plot the result
% plotZono(cZono);

% compare with ground-truth
V_ = [-1 0 -2;3 0 -2];

if ~compareMatrices(V,V_)
    res = false;
end



% TEST 3 ------------------------------------------------------------------

% construct zonotope
Z = [0 3 0 1 -2;0 0 2 1 1];
A = [0 0 0 1];
b = 0.5;
cZono = conZonotope(Z,A,b);

% calculate vertices
V = vertices(cZono);

% plot the result
% plotZono(cZono);

% compare with ground-truth
V_ = [-4 -4 -2 4 4 2;-3 1 3 3 -1 -3] + [-1;0.5]*ones(1,6);

if ~compareMatrices(V,V_)
    res = false;
end


% TEST 4 ------------------------------------------------------------------

% construct zonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZono = conZonotope(Z,A,b);

% calculate vertices
V = vertices(cZono);

% plot the result
% plotZono(cZono);

% compare with ground-truth
V_ = [1 1 3 3;3 -1 2 -2];

if ~compareMatrices(V,V_)
    res = false;
end

%------------- END OF CODE --------------
