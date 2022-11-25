function matZ = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or a 
% matrix zonotope with a matrix zonotope
%
% Syntax:  
%    matZ = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical matrix or matrix zonotope
%    factor2 - numerical matrix or matrix zonotope
%
% Outputs:
%    matZ - matrix zonotope
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
% Last update:  05-August-2010
% Last revision:---

%------------- BEGIN CODE --------------

%factor1 is a numeric matrix
if isnumeric(factor1)
    %initialize factor
    matrix=factor1;
    %initialize matrix zonotope
    matZ=factor2;
    %compute center
    matZ.center=matrix*matZ.center;
    for i=1:matZ.gens
        matZ.generator{i}=matrix*matZ.generator{i};
    end
    
%factor2 is a numeric matrix
elseif isnumeric(factor2)
    %initialize factor
    matrix=factor2;
    %initialize matrix zonotope
    matZ=factor1;
    %compute center
    matZ.center=matZ.center*matrix;
    for i=1:matZ.gens
        matZ.generator{i}=matZ.generator{i}*matrix;
    end
    
%both factors are zonotope matrices
else
    %initialize matrix zonotope
    matZ1=factor1;
    %initialize matrix zonotope
    matZ2=factor2;
    %initialize matrix zonotope
    matZ=matZonotope();
    %compute center
    matZ.center=matZ1.center*matZ2.center;
    %compute generators
    %center1 with generators2
    for i=1:matZ2.gens
        matZ.generator{end+1}=matZ1.center*matZ2.generator{i};
    end
    %generator1 with center2
    for i=1:matZ1.gens
        matZ.generator{end+1}=matZ1.generator{i}*matZ2.center;
    end
    %generators1 with generators2
    for j=1:matZ1.gens
        for i=1:matZ2.gens
            matZ.generator{end+1}=matZ1.generator{j}*matZ2.generator{i};
        end
    end
    %update number of generators and dimension
    matZ.dim=matZ.dim;
    matZ.gens=(matZ1.gens+1)*(matZ2.gens+1)-1;
end


%------------- END OF CODE --------------