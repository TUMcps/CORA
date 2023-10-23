function res = test_conZonotope_isemptyobject
% test_conZonotope_isemptyobject - unit test function of isemptyobject
%
% Syntax:
%    res = test_conZonotope_isemptyobject
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
% Written:       03-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate constrained zonotopes
cZ1 = conZonotope();

Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1];
b = 1;
cZ2 = conZonotope(Z,A,b);

% check results
res = isemptyobject(cZ1) && ~isemptyobject(cZ2);

% ------------------------------ END OF CODE ------------------------------
