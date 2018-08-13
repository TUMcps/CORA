function [Gnew]=additionalGenerators(dim,sections)
% additionalGenerators - computes evenly distributed directions in all
% dimensions
%
% Syntax:  
%    [G]=additionalGenerators(dim,sections)
%
% Inputs:
%    dim - dimension
%    sections - number of sections in each dimension
%
% Outputs:
%    G - matrix of generators/directions
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      19-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%first unit vector
Gnew = [1;zeros(dim-1,1)];

%increment
inc = pi/(sections+1);

%loop over all dimensions
for iDim = 1 : dim-1
    %set up Gnew
    G = Gnew;
    Gnew = [];
    for iAngle = 1 : sections
        %determine angle
        angle = iAngle*inc;
        
        %construct rotational matrix
        RotMat = eye(dim);
        RotMat(iDim, iDim) = cos(angle);
        RotMat(iDim, iDim+1) = -sin(angle);
        RotMat(iDim+1, iDim) = sin(angle);
        RotMat(iDim+1, iDim+1) = cos(angle);
        
        %compute new generators
        for iGen = 1 : length(G(1,:))
            Gnew(:,end+1) = RotMat*G(:,iGen);
        end
    end
end



%------------- END OF CODE --------------