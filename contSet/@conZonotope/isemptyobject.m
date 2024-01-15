function res = isemptyobject(cZ)
% isemptyobject - checks whether a constrained zonotope contains any
%    information at all; consequently, the set is interpreted as the empty
%    set 
%
% Syntax:
%    res = isemptyobject(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    res - true/false
%
% Example: 
%    cZ = conZonotope([2;1],[1;2],1,0);
%    isemptyobject(cZ); % false
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

res = false(size(cZ));
% loop over class-arrays
for i=1:size(cZ,1)
    for j=1:size(cZ,2)
        res(i,j) = aux_checkIfEmpty(cZ(i,j));
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkIfEmpty(cZ)

    res = isnumeric(cZ.c) && isempty(cZ.c) ...
        &&isnumeric(cZ.G) && isempty(cZ.G) ...
        && isnumeric(cZ.A) && isempty(cZ.A) ...
        && isnumeric(cZ.b) && isempty(cZ.b) ...
        && isnumeric(cZ.ksi) && isempty(cZ.ksi) ...
        && isnumeric(cZ.R) && isempty(cZ.R);

end

% ------------------------------ END OF CODE ------------------------------
