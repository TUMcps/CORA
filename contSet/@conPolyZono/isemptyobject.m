function res = isemptyobject(cPZ)
% isemptyobject - checks whether a constrained polynomial zonotope contains
%    any information at all; consequently, the set is interpreted as the
%    empty set 
%
% Syntax:
%    res = isemptyobject(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object
%
% Outputs:
%    res - true/false
%
% Example: 
%    c = [0;0];
%    G = [1 0 1 -1; 0 1 1 1];
%    E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
%    A = [1 -0.5 0.5];
%    b = 0.5;
%    EC = [0 1 2; 1 0 0; 0 1 0];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%    res = isemptyobject(cPZ); % false
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

res = false(size(cPZ));
% loop over class-arrays
for i=1:size(cPZ,1)
    for j=1:size(cPZ,2)
        res(i,j) = aux_checkIfEmpty(cPZ(i,j));
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkIfEmpty(cPZ)

    res = isnumeric(cPZ.c) && isempty(cPZ.c) ...
        && isnumeric(cPZ.G) && isempty(cPZ.G) ...
        && isnumeric(cPZ.E) && isempty(cPZ.E) ...
        && isnumeric(cPZ.A) && isempty(cPZ.A) ...
        && isnumeric(cPZ.b) && isempty(cPZ.b) ...
        && isnumeric(cPZ.EC) && isempty(cPZ.EC) ...
        && isnumeric(cPZ.GI) && isempty(cPZ.GI) ...
        && isnumeric(cPZ.id) && isempty(cPZ.id);

end

% ------------------------------ END OF CODE ------------------------------
