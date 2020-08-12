function intMatPower = mpower(intMat,exponent)
% mpower - Overloaded '^' operator for the power of an interval matrix
%
% Syntax:  
%    intMatPower = mpower(intMat,exponent)
%
% Inputs:
%    intMat - interval matrix 
%    exponent - exponent
%
% Outputs:
%    intMatPower - interval matrix 
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
% Last update:  05-August-2010
% Last revision:---

%------------- BEGIN CODE --------------

%factor1 is a numeric matrix
if exponent>=0
    if exponent==0
        %return identity matrix
        intMatPower=intMat;
        intMatPower.int=intMat.int^0;
    elseif exponent==1
        %do nothing
        intMatPower=intMat;
    else
        intMatPower=intMat*intMat;
        for i=3:exponent
        %multiply matrix zonotope with itself
            intMatPower=intMatPower*intMat;
        end
    end
    
else
    error('no negative powers supported')
end

%------------- END OF CODE --------------