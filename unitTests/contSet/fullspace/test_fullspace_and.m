function res = test_fullspace_and
% test_fullspace_and - unit test function of and
%
% Syntax:
%    res = test_fullspace_and
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
% Written:       25-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init fullspace
n = 2;
fs = fullspace(n);

% intersection with itself
assert(isequal(fs & fs,fs));

% intersection with zonotope
Z = zonotope(zeros(n,1),eye(n));
assert(isequal(fs & Z,Z));

% intersection with emptySet
O = emptySet(n);
assert(isequal(fs & O,O));

% unbounded interval
I = interval([-Inf;2],[9;Inf]);
assert(isequal(fs & I,I));

% numeric vector
p = [1;-1];
assert(compareMatrices(fs & p,p));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
