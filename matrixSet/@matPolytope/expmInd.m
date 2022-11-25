function eP = expmInd(matP,maxOrder)
% expmInd - operator for the exponential matrix of a 
% matrix polytope, evaluated independently
%
% Syntax:  
%    eP = expmInd(matP,maxOrder)
%
% Inputs:
%    matP - matrix polytope
%    maxOrder - maximum Taylor series order until remainder is computed
%
% Outputs:
%    eP - matrix polytope exponential
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

%compute finite Taylor series
%initialize matrix zonotope
eP=matP^0;
%initialize power
ePpow=matP^0;

%compute finite Taylor sum
for i=1:maxOrder
    ePpow = ePpow*matP;
    eP = eP + ePpow*(1/factorial(i));
end

%compute remainder value
%create over-approximating interval matrix
intMat = intervalMatrix(matP);
%compute remainder value
E = exponentialRemainder(intMat,maxOrder);

%convert remainder and add it to the Taylor series
eP = eP + matPolytope(E);

%------------- END OF CODE --------------