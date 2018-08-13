function Z = enlarge(Z,factorVec)
% enlarge - enlarges the generators of a zonotope by a vector of factors
%
% Syntax:  
%    Z = enlarge(Z,factorVec)
%
% Inputs:
%    Z - zonotope object
%    factorVec - vector of factors for the enlargement of each dimension
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      20-November-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

Z.Z(:,2:end)=diag(factorVec)*Z.Z(:,2:end);

%------------- END OF CODE --------------