function S_out = decompose(S,blocks)
% decompose - block decomposition of a set into a cell array of projected
%    sets (no cell if only one block)
%
% Syntax:
%    S_out = decompose(S,blocks)
%
% Inputs:
%    S - contSet object
%    blocks - bx2 array with b blocks
%
% Outputs:
%    S_out - cell array of contSet objects
%
% Example:
%    S = zonotope([1;-1;0;0;1],eye(5));
%    blocks = [1 2; 3 5];
%    S_out = decompose(S,blocks);
%
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reach_decomp

% Authors:       Mark Wetzlinger
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% number of blocks
numBlocks = size(blocks,1);

% no projection if a single block
if numBlocks == 1
    S_out = S.copy();
    return
end

% use project function
S_out = arrayfun(@(i) project(S,blocks(i,1):blocks(i,2)),...
                 1:numBlocks,'UniformOutput',false)';

% ------------------------------ END OF CODE ------------------------------
