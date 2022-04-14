function matPpower = mpower(matP,exponent)
% mpower - Overloaded '^' operator for the power of matrix polytope 
%
% Syntax:  
%    matPpower = mpower(matP,exponent)
%
% Inputs:
%    matP - matrix polytope
%    exponent - exponent
%
% Outputs:
%    matP - matrix polytope
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      21-June-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%factor1 is a numeric matrix
if exponent>=0
    if exponent==0
        %return identity matrix
        matPpower=matPolytope();
        matPpower.dim=matP.dim;
        matPpower.verts=1;
        matPpower.vertex{1}=eye(matP.dim);       
    elseif exponent==1
        %do nothing
        matPpower=matP;
    else
        matPpower=matP*matP;
        for i=3:exponent
        %multiply matrix zonotope with itself
            matPpower=matPpower*matP;
        end
    end
else
    matPpower=[];
    error('no negative powers supported')
end

%------------- END OF CODE --------------