function res = isemptyobject(pZ)
% isemptyobject - checks whether a polynomial zonotope contains any
%    information at all; consequently, the set is interpreted as the empty
%    set
%
% Syntax:
%    res = isemptyobject(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    res - true/false
%
% Example: 
%    pZ = polyZonotope([1;-1],[0;1]);
%    isemptyobject(pZ); % false
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

res = false(size(pZ));
% loop over class-arrays
for i=1:size(pZ,1)
    for j=1:size(pZ,2)
        res(i,j) = aux_checkIfEmpty(pZ(i,j));
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkIfEmpty(pZ)

    res = isnumeric(pZ.c) && isempty(pZ.c) ...
        && isnumeric(pZ.G) && isempty(pZ.G) ...
        && isnumeric(pZ.GI) && isempty(pZ.GI) ...
        && isnumeric(pZ.E) && isempty(pZ.E) ...
        && isnumeric(pZ.id) && isempty(pZ.id);

end

% ------------------------------ END OF CODE ------------------------------
