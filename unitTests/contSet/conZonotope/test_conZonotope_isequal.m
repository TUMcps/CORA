function res = test_conZonotope_isequal
% test_conZonotope_isequal - unit test function for set equality check
%
% Syntax:
%    res = test_conZonotope_isequal
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

% Authors:       Mark Wetzlinger
% Written:       19-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% constrained zonotope 1
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZ1 = conZonotope(Z,A,b);

% check equality with itself
if ~isequal(cZ1,cZ1)
    res = false;
end

% shift by a small offset
if isequal(cZ1,cZ1+1e-8*ones(2,1))
    res = false;
end

% compensate for small offset by tolerance
if ~isequal(cZ1,cZ1+1e-8*ones(2,1),1e-6)
    res = false;
end

% shift constraints by a small offset
b_ = 1+eps;
cZ1_ = conZonotope(Z,A,b_);

if ~isequal(cZ1,cZ1_)
    res = false;
end

% completely different constrained zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1];
b = 1;
cZ2 = conZonotope(Z,A,b);

% should not be equal
if isequal(cZ1,cZ2)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
