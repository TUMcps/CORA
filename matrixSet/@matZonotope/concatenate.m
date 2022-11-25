function matZ1 = concatenate(matZ1,matZ2)
% concatenate - concatenates the center and all generators of the second
% matrix zonotope to the first one
%
% Syntax:  
%    matZ = concatenate(matZ1,matZ2)
%
% Inputs:
%    matZ1 - zonotope matrix object
%    matZ2 - zonotope matrix object
%
% Outputs:
%    matZ - zonotope matrix object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      06-September-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%add center
matZ1.generator{end+1} = matZ2.center;

%add generators
for i = 1:length(matZ2.generator)
    matZ1.generator{end+1} = matZ2.generator{i};
end

%update number of generators
matZ1.gens = length(matZ1.generator);


%------------- END OF CODE --------------