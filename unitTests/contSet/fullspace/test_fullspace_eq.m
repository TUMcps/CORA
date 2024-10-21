function res = test_fullspace_eq
% test_fullspace_eq - unit test function of box
%
% Syntax:
%    res = test_fullspace_eq
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
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init fullspace
n = 2;
fs = fullspace(n);

% compare with itself
assert(fs == fs);

% init zonotope
Z = zonotope(zeros(n,1),eye(n));
assert(~(fs == Z));

% init interval
I = interval(-Inf(n,1),Inf(n,1));
assert(fs == I);

% init capsule
C = capsule([1;1],[-1;1],0.5);
assert(~(fs == C));

% numerical vector
p = [2;1];
assert(~(fs == p));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
