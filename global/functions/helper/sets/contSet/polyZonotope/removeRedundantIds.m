function [E_unique,id_unique] = removeRedundantIds(E,id)
% removeRedundantIds - finds and removes redundant elements in the
% ID-vectors and sums the corresponding rows of the exponent matrix.
%
% Syntax:
%    [id,E1,E2] = removeRedundantIds(E,id)
%
% Inputs:
%    id - ID-vector of the polynomial zonotope
%    E  - exponent matrix of the polynomial zonotope
%
% Outputs:
%    id - merged ID-vector
%    E1 - adapted exponent matrix of the first polynomial zonotope
%    E2 - adapted exponent matrix of the second polynomial zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Lukas Schäfer
% Written:       07-February-2025 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % find unique ids (+ indices in ID-vector)
    id_unique = unique(id, 'stable');

    % check if exponent matrix is empty
    if isempty(E)
        E_unique = [];
        return;
    end
    
    % initialize unique exponent matrix
    E_unique = zeros(length(id_unique),size(E,2));
    % iterate over unique ids to sum corresponding exponents
    for ii = 1:length(id_unique)
        E_unique(ii,:) = sum(E(id == id_unique(ii),:),1);
    end

end

% ------------------------------ END OF CODE ------------------------------
