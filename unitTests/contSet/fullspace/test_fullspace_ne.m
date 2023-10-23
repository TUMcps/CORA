function res = test_fullspace_ne
% test_fullspace_ne - unit test function of ne
%
% Syntax:
%    res = test_fullspace_ne
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
res = ~(fs ~= fs);

% init zonotope
Z = zonotope(zeros(n,1),eye(n));
res(end+1,1) = fs ~= Z;

% init interval
I = interval(-Inf(n,1),Inf(n,1));
res(end+1,1) = ~(fs ~= I);

% init capsule
C = capsule([1;1],[-1;1],0.5);
res(end+1,1) = fs ~= C;

% numerical vector
p = [2;1];
res(end+1,1) = fs ~= p;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
