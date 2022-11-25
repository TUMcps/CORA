function Z = zonotope(matZ)
% zonotope - Converts a matrix zonotope into a zonotope 
%
% Syntax:  
%    Z = zonotope(matZ)
%
% Inputs:
%    matZ - matrix zonotope
%
% Outputs:
%    Z - zonotope
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      18-June-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%concert center
center=mat2vec(matZ.center);

%convert generators
for i=1:matZ.gens
    generatorMatrix(:,i) = mat2vec(matZ.generator{i});
end
    
%instantiate zonotope
Z=zonotope([center,generatorMatrix]);

%------------- END OF CODE --------------