function S_out = block_mtimes(M,S)
% block_mtimes - evaluation of a linear map using block decomposition
%
% Syntax:
%    S_out = block_mtimes(M,S)
%
% Inputs:
%    M - numeric matrix, interval object
%    S - cell array of sets, contSet object, numeric vector
%
% Outputs:
%    S_out - sets after mapping
%
% Example:
%    M = eye(5);
%    S = {zonotope(ones(2,1)),zonotope(-ones(3,1))};
%    S = block_mtimes(M,S);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/priv_reach_decomp

% Authors:       Mark Wetzlinger
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{M,'att',{'numeric','interval'}};...
                {S,'att',{'cell','contSet','numeric'}}});

% quick exit for contSet and numeric
if isa(S,'contSet') || isnumeric(S)
    S_out = M*S;
    return
end

% read out block sizes
numBlocks = numel(S);
if isnumeric(S{1})
    dims = cellfun(@(S_) length(S_),S,'UniformOutput',true);
else
    dims = cellfun(@(S_) S_.dim,S,'UniformOutput',true);
end
blocks_end = reshape(cumsum(dims),[],1);
blocks = [[1; blocks_end(1:end-1)+1], blocks_end];

% init new cell array for result
S_out = cell(numel(S),1);

% formula: \forall i: S_i = \sum_j M_ij * S_j
for i=1:numBlocks
    S_out{i} = zeros(dims(i),1);
    for j=1:numBlocks
        M_block = M(blocks(i,1):blocks(i,2),blocks(j,1):blocks(j,2));
        S_out{i} = S_out{i} + M_block * S{j};
    end
end

% ------------------------------ END OF CODE ------------------------------
