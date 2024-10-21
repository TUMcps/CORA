function res = test_block_mtimes()
% test_block_mtimes - unit test function for block-wise multiplication of a
%    matrix and a set
%
% Syntax:
%    res = test_block_mtimes()
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
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-12;

% init matrix and sets
M = [1  2  3 -2 1;
     0 -1  2  0 0;
     1  1  4 -1 4;
     3  0 -2  2 0;
    -1  1  2  0 5];
S = interval([-1;2;0;0;2],[4;2;1;0;5]);
S_cell{1,1} = interval([-1;2],[4;2]);
S_cell{2,1} = interval([0;0;2],[1;0;5]);

% single block
S_mtimes = block_mtimes(M,S);
S_true = M*S;
assert(isequal(S_mtimes,S_true,tol));

% mulitple blocks
S_mtimes = block_mtimes(M,S_cell);
assert(iscell(S_mtimes));
S_true = cell(2,1);
S_true{1} = M(1:2,1:2) * S_cell{1} + M(1:2,3:5) * S_cell{2};
S_true{2} = M(3:5,1:2) * S_cell{1} + M(3:5,3:5) * S_cell{2};
assert(isequal(S_mtimes{1},S_true{1},tol));
assert(isequal(S_mtimes{2},S_true{2},tol));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
