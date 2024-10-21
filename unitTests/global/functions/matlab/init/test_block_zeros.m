function res = test_block_zeros()
% test_block_zeros - unit test function for instantiation of all-zero
%    vectors in block-specific sizes
%
% Syntax:
%    res = test_block_zeros()
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

% single block
blocks = [1 5];
S = block_zeros(blocks);
assert(compareMatrices(S,zeros(5,1),0,"equal",true));

% multiple blocks
blocks = [1 2; 3 5; 6 10];
S = block_zeros(blocks);
assert(iscell(S));
assert(compareMatrices(S{1},zeros(2,1),0,"equal",true));
assert(compareMatrices(S{2},zeros(3,1),0,"equal",true));
assert(compareMatrices(S{3},zeros(5,1),0,"equal",true));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
