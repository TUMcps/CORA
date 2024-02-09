function simRes = uminus(simRes)
% uminus - Overloads the unary '-' operator
%
% Syntax:
%    simRes = -simRes
%    simRes = uminus(simRes)
%
% Inputs:
%    simRes - simResult object
%
% Outputs:
%    simRes - transformed simResult object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       02-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

simRes = -1 * simRes;

end

% ------------------------------ END OF CODE ------------------------------
