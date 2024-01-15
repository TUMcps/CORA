function res = isemptyobject(probZ)
% isemptyobject - checks whether a probabilistic zonotope contains any
%    information at all; consequently, the set is interpreted as the empty
%    set 
%
% Syntax:
%    res = isemptyobject(probZ)
%
% Inputs:
%    probZ - probZonotope object
%
% Outputs:
%    res - true/false
%
% Example: 
%    probZ = probZonotope([10 1 -2; 0 1 1],[0.6 1.2; 0.6 -1.2]);
%    isemptyobject(probZ); % false
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

res = false(size(probZ));
% loop over class-arrays
for i=1:size(probZ,1)
    for j=1:size(probZ,2)
        res(i,j) = aux_checkIfEmpty(probZ(i,j));
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkIfEmpty(probZ)

    res = isnumeric(probZ.Z) && isempty(probZ.Z) ...
        && isnumeric(probZ.g) && isempty(probZ.g) ...
        && isnumeric(probZ.cov) && isempty(probZ.cov) ...
        && islogical(probZ.gauss) && ~probZ.gauss ...
        && isnumeric(probZ.gamma) && isscalar(probZ.gamma) && probZ.gamma == 2;

end

% ------------------------------ END OF CODE ------------------------------
