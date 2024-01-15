function res = isemptyobject(Z)
% isemptyobject - checks whether a zonotope contains any information at
%    all; consequently, the set is interpreted as the empty set 
%
% Syntax:
%    res = isemptyobject(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    res - true/false
%
% Example: 
%    Z = zonotope([2;1]);
%    isemptyobject(Z); % false
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

res = false(size(Z));
% loop over class-arrays
for i=1:size(Z,1)
    for j=1:size(Z,2)
        res(i,j) = aux_checkIfEmpty(Z(i,j));
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkIfEmpty(Z)

    res = isnumeric(Z.c) && isempty(Z.c) ...
        && isnumeric(Z.G) && isempty(Z.G) ...
        && isnumeric(Z.halfspace) && isempty(Z.halfspace);

end

% ------------------------------ END OF CODE ------------------------------
