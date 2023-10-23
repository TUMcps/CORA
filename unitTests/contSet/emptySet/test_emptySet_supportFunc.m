function res = test_emptySet_supportFunc
% test_emptySet_supportFunc - unit test function of and
%
% Syntax:
%    res = test_emptySet_supportFunc
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

% init empty set
n = 2;
O = emptySet(n);

% direction
dir = [1;1];

% upper
[val,x] = supportFunc(O,dir);
res = val == -Inf && all(size(x) == [n,0]);

% lower
[val,x] = supportFunc(O,dir,'lower');
res(end+1,1) = val == Inf && all(size(x) == [n,0]);

% range
val = supportFunc(O,dir,'range');
res(end+1,1) = val == interval(-Inf,Inf);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
