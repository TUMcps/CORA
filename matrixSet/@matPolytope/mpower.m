function matPpower = mpower(matP,exponent)
% mpower - Overloaded '^' operator for the power of matrix polytope 
%
% Syntax:
%    matPpower = mpower(matP,exponent)
%
% Inputs:
%    matP - matPolytope object
%    exponent - exponent
%
% Outputs:
%    matPpower - ?
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       21-June-2010 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{matP,'att','matPolytope'}, ...
                {exponent,'att','numeric',{'integer','nonnegative','scalar'}}});

%factor1 is a numeric matrix
if exponent==0
    %return identity matrix
    matPpower=matPolytope();
    matPpower.verts=1;
    matPpower.vertex{1}=eye(dim(matP));

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

% ------------------------------ END OF CODE ------------------------------
