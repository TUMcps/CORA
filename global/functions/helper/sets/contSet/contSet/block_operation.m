function S_out = block_operation(op,varargin)
% block_operation - evaluates a set operation over a cell array of sets,
%    where each pair must be of equal dimensions; the primary purpose
%    of this helper function is to ease the reading of the computation of
%    auxiliary terms for the reachable set computation of linear systems
%
% Syntax:
%    S_out = block_operation(op,S1)
%    S_out = block_operation(op,S1,S2,...)
%
% Inputs:
%    op - function handle to a binary set operation
%    S1 - contSet object, cell array of sets
%    S2 - contSet object, cell array of sets
%
% Outputs:
%    S_out - (cell array of) result(s)
%
% Example:
%    S1 = {zonotope(0),zonotope([1;-1])};
%    S2 = {zonotope(1),zonotope([4;2])};
%    S_out = block_operation(@plus,S1,S2);
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

% check if all sets are contSet or numeric types -> evaluate as normally
if all(cellfun(@(S) isa(S,'contSet') || isnumeric(S),...
               varargin,'UniformOutput',true))
    S_out = op(varargin{:});
    return
end

% otherwise, all input arguments must be sets
if all(cellfun(@(S) iscell(S),varargin,'UniformOutput',true))
    % take number of elements of first group of sets
    numOperands = numel(varargin);
    numBlocks = numel(varargin{1});
    % evaluate operation for each pair of sets... reshape
    S_tuplets = arrayfun(@(i) arrayfun(@(j) varargin{j}{i},1:numOperands,'UniformOutput',false), 1:numBlocks, 'UniformOutput', false);
    S_out = cellfun(@(S) op(S{:}),S_tuplets,'UniformOutput',false);
    return
end

% ------------------------------ END OF CODE ------------------------------
