function P = polytope(matP)
% polytope - Converts a matrix polytope to a polytope 
%
% Syntax:
%    P = polytope(matP)
%
% Inputs:
%    matP - matrix polytope
%
% Outputs:
%    P - polytope
%
% Example: 
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       21-June-2010
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert vertices
vecV = reshape(matP.V,[],matP.numverts);

% instantiate polytope 
P=polytope(vecV);

% ------------------------------ END OF CODE ------------------------------
