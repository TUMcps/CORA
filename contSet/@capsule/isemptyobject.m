function res = isemptyobject(C)
% isemptyobject - checks whether a capsule contains any information at all;
%    consequently, the set is interpreted as the empty set 
%
% Syntax:
%    res = isemptyobject(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    res - true/false
%
% Example: 
%    C = capsule([1;-1],[0;1],1);
%    isemptyobject(C); % false
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

res = false(size(C));
% loop over class-arrays
for i=1:size(C,1)
    for j=1:size(C,2)
        res(i,j) = aux_checkIfEmpty(C(i,j));
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkIfEmpty(C)

    res = isnumeric(C.c) && isempty(C.c) ...
        && isnumeric(C.g) && isempty(C.g) ...
        && isnumeric(C.r) && isempty(C.r);

end

% ------------------------------ END OF CODE ------------------------------
