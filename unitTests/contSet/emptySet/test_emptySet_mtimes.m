function res = test_emptySet_mtimes
% test_emptySet_mtimes - unit test function of mtimes
%
% Syntax:
%    res = test_emptySet_mtimes
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

% multiplication with scalar
s = 2;
O_ = s*O;
res = isequal(O,O_);
O_ = O*s;
res(end+1,1) = isequal(O,O_);

% multiplication with a square matrix
M = [2 1; -1 3];
O_ = M*O;
res(end+1,1) = isequal(O,O_);

% multiplication as projection onto a subspace
M = [2 1];
O_ = M*O;
res(end+1,1) = isequal(emptySet(1),O_);

% multiplication as projection to a higher-dimensional space
M = [2 1; 1 3; -1 0];
O_ = M*O;
res(end+1,1) = isequal(emptySet(3),O_);

% multiplication with an interval matrix
intMat = intervalMatrix([2 1],[1 1]);
O_ = intMat*O;
res(end+1,1) = isequal(emptySet(1),O_);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
