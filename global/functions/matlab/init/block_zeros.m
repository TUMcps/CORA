function S_out = block_zeros(blocks)
% block_zeros - initializes a cell array for a number of blocks with
%    all-zero column vectors of block dimension; if a single block is
%    given, we instead return a vector instead of a cell
%
% Syntax:
%    S_out = block_zeros(blocks)
%
% Inputs:
%    blocks - bx2 array with b blocks
%
% Outputs:
%    S_out - (cell array of) all-zero vector(s) of block dimension
%
% Example:
%    S_out = block_zeros([1 2; 3 5]);
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

numBlocks = size(blocks,1);
if numBlocks == 1
    S_out = zeros(blocks(end,end),1);
else
    S_out = arrayfun(@(i) zeros(blocks(i,2)-blocks(i,1)+1,1),1:numBlocks,...
        'UniformOutput',false)';
end

% ------------------------------ END OF CODE ------------------------------
