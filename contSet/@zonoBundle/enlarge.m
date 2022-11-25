function Zbundle = enlarge(Zbundle,factorVec)
% enlarge - enlarges the generators of a zonotope bundle by a vector of factors
%
% Syntax:  
%    Zbundle = enlarge(Zbundle,factorVec)
%
% Inputs:
%    Zbundle - zonotope bundle
%    factorVec - vector of factors for the enlargement of each dimension
%
% Outputs:
%    Zbundle - zonotope bundle
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

for i=1:Zbundle.parallelSets
    Zbundle.Z{i}=enlarge(Zbundle.Z{i},factorVec);
end

%------------- END OF CODE --------------