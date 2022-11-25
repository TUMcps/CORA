function matZpower = mpower(matZ,exponent)
% mpower - Overloaded '^' operator for the power of matrix zonotope 
%
% Syntax:  
%    matZ = mpower(matZ,exponent)
%
% Inputs:
%    matZ - matrix zonotope
%    exponent - exponent
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
if exponent>=0
    if exponent==0
        %return identity matrix
        matZpower=matZ;
        matZpower.center=eye(matZ.dim);
        matZpower.generator=[];
        matZpower.gens=0;
    elseif exponent==1
        %do nothing
        matZpower=matZ;
    else
        matZpower=matZ*matZ;
        for i=3:exponent
        %multiply matrix zonotope with itself
            matZpower=matZpower*matZ;
        end
    end
else
    matZpower=[];
    disp('no negative powers supported')
end

%------------- END OF CODE --------------