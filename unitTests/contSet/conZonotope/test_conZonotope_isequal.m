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
assert(isequal(cZ1,cZ1))


% shift by a small offset
assert(~isequal(cZ1,cZ1+1e-8*ones(2,1)))

% compensate for small offset by tolerance
assert(isequal(cZ1,cZ1+1e-8*ones(2,1),1e-6))

% shift constraints by a small offset
b_ = 1+eps;
cZ1_ = conZonotope(Z,A,b_);

assert(isequal(cZ1,cZ1_))

% completely different constrained zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1];
b = 1;
cZ2 = conZonotope(Z,A,b);

% should not be equal
assert(~isequal(cZ1,cZ2))

end

% ------------------------------ END OF CODE ------------------------------
