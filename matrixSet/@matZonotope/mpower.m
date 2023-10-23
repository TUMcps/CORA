function matZpower = mpower(matZ,exponent)
% mpower - Overloaded '^' operator for the power of matrix zonotope 
%
% Syntax:
%    matZ = mpower(matZ,exponent)
%
% Inputs:
%    matZ - matZonotope object
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

% Authors:       Matthias Althoff
% Written:       18-June-2010 
% Last update:   05-August-2010
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{matZ,'att','matZonotope'}, ...
                {exponent,'att','numeric',{'scalar','integer','nonnegative'}}});

%factor1 is a numeric matrix
if exponent==0
    %return identity matrix
    matZpower=matZ;
    matZpower.center=eye(dim(matZ));
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

% ------------------------------ END OF CODE ------------------------------
