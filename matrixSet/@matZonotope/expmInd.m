function eZ = expmInd(matZ,maxOrder)
% expmInd - operator for the exponential matrix of a 
% matrix zonotope, evaluated independently
%
% Syntax:  
%    eZ = expmInd(matZ,maxOrder)
%
% Inputs:
%    matZ - matrix zonotope
%    maxOrder - maximum Taylor series order until remainder is computed
%
% Outputs:
%    eZ - matrix zonotope exponential
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

%compute finite Taylor series
%initialize matrix zonotope
eZ=matZ^0;
%initialize power
eZpow=matZ^0;

%compute finite Taylor sum
for i=1:maxOrder
    eZpow = eZpow*matZ;
    eZ = eZ + eZpow*(1/factorial(i));
end

%compute remainder value
%create over-approximating interval matrix
intMat = intervalMatrix(matZ);
%compute remainder value
E = exponentialRemainder(intMat,maxOrder);

%convert remainder and add it to the Taylor series
eZ = eZ + matZonotope(E);

%------------- END OF CODE --------------