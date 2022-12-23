function res = test_conZonotope_conZonotope
% test_conZonotope_conZonotope - unit test function of conZonotope (constructor)
%
% Syntax:  
%    res = test_conZonotope_conZonotope
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

% Author:       Mark Wetzlinger
% Written:      19-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty conZonotope
cZ = conZonotope();
res = true;
if ~isempty(cZ)
    res = false;
end

% init simple constrained zonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZ = conZonotope(Z,A,b);

if ~all(all(withinTol(cZ.Z,Z))) || ~all(all(withinTol(cZ.A,A))) ...
        || ~all(withinTol(cZ.b,b))
    res = false;
end

%------------- END OF CODE --------------