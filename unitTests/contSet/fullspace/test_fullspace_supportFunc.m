function res = test_fullspace_supportFunc
% test_fullspace_supportFunc - unit test function of supportFunc
%
% Syntax:
%    res = test_fullspace_supportFunc
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

% direction
dir = [1;0.5];

% compute support function
[val,x] = supportFunc(fs,dir);
assert(val == Inf && all(x == [Inf;Inf]));
[val,x] = supportFunc(fs,dir,'lower');
assert(val == -Inf && all(x == [-Inf;-Inf]));
[val,x] = supportFunc(fs,dir,'range');
assert(val == interval(-Inf,Inf) && all(all(x == [-Inf Inf;-Inf Inf])));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
