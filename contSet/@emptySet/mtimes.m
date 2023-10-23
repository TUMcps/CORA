function O = mtimes(factor1,factor2)
% mtimes - overloaded '*' operator for the linear map of an empty set
%
% Syntax:
%    O = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - emptySet object, numerical scalar/matrix
%    factor2 - emptySet object, numerical scalar/matrix
%
% Outputs:
%    O - emptySet object
%
% Example: 
%    O = emptySet(2);
%    M = [2 1; -1 3];
%    M*O
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   05-April-2023 (MW, bug fix)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% find the fullspace object
[O,M] = findClassArg(factor1,factor2,'emptySet');

if isscalar(M)
    % keep O as is...
    
elseif isnumeric(M)
    % projection to subspace or higher-dimensional space
    O.dimension = size(M,1);

elseif isa(M,'intervalMatrix') || isa(M,'matZonotope') || isa(M,'matPolytope')
    % projection to subspace or higher-dimensional space 
    O.dimension = dim(M,1);

end

% ------------------------------ END OF CODE ------------------------------
