function res = isemptyobject(E)
% isemptyobject - checks whether an ellipsoid contains any information at
%    all; consequently, the set is interpreted as the empty set 
%
% Syntax:
%    res = isemptyobject(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    res - true/false
%
% Example: 
%    E = ellipsoid([1,0;0,2],[0;1]);
%    isemptyobject(E); % false
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

res = false(size(E));
% loop over class-arrays
for i=1:size(E,1)
    for j=1:size(E,2)
        res(i,j) = aux_checkIfEmpty(E(i,j));
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkIfEmpty(E)

    res = isnumeric(E.Q) && isempty(E.Q) ...
        && isnumeric(E.q) && isempty(E.q) ...
        && isnumeric(E.TOL) && ...
        ( (isscalar(E.TOL) && E.TOL == 1e-6) || isempty(E.TOL) );

end

% ------------------------------ END OF CODE ------------------------------
