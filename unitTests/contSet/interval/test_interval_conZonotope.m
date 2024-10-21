function res = test_interval_conZonotope
% test_interval_conZonotope - unit test function of conversion to
%    conZonotopes
%
% Syntax:
%    res = test_interval_conZonotope
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
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init interval
I = interval([-2;-1],[4;7]);

% convert to constrained zonotope
cZ = conZonotope(I);

% check equality of vertices
V_I = vertices(I);
V_cZ = vertices(cZ);
assert(compareMatrices(V_I,V_cZ));
% check set equality
assert(isequal(cZ,I));


% degenerate interval
I = interval([-2;0;5],[-1;0;7]);

% convert to constrained zonotope
cZ = conZonotope(I);

% check equality of vertices
V_I = vertices(I);
V_cZ = vertices(cZ);
assert(compareMatrices(V_I,V_cZ));
% check set equality
assert(isequal(cZ,I));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
