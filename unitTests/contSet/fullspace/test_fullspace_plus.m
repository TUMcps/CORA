function res = test_fullspace_plus
% test_fullspace_plus - unit test function of plus
%
% Syntax:
%    res = test_fullspace_plus
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

% init vector
p = [2;1];
fs_ = fs + p;
assert(isequal(fs,fs_));

% empty vector
p = double.empty(n,0);
fs_ = fs + p;
assert(isequal(fs_,emptySet(n)));

% init zonotope
Z = zonotope(zeros(n,1),eye(n));
fs_ = fs + Z;
assert(isequal(fs,fs_));

% init interval
I = interval([-2;1],[Inf;3]);
fs_ = fs + I;
assert(isequal(fs,fs_));

% init emptySet
O = emptySet(n);
fs_ = fs + O;
assert(isequal(fs_,O));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
