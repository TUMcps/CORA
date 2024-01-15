function res = isemptyobject(zB)
% isemptyobject - checks whether a zonotope bundle contains any information
%    at all; consequently, the set is interpreted as the empty set 
%
% Syntax:
%    res = isemptyobject(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    res - true/false
%
% Example: 
%    zB = zonoBundle({zonotope([1;0],[1 0; 0 1]),zonotope([-1;2],[1 1; -1 0])});
%    isemptyobject(zB); % false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       24-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = false(size(zB));
% loop over class-arrays
for i=1:size(zB,1)
    for j=1:size(zB,2)
        res(i,j) = aux_checkIfEmpty(zB(i,j));
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkIfEmpty(zB)

    res = (iscell(zB.Z) && isempty(zB.Z) ...
        && isnumeric(zB.parallelSets) && isscalar(zB.parallelSets) && zB.parallelSets == 0) ...
        || all(cellfun(@isemptyobject,zB.Z,'UniformOutput',true));

end

% ------------------------------ END OF CODE ------------------------------
