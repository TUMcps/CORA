function scaling = scalingSquaringHornerTaylorSeries(intMat,maxOrder,potentiation)
% returns the approximation of e^intMat using different algorithms 
%    with maxOrder iterations. It is used as an fassade to access the
%    algorithms in the private directory
%
% Syntax:  
%    scaling = scalingSquaringHornerTaylorSeries(intMat,maxOrder,potentiation)
%
% Inputs:
%    intMat - intervalMatrix (nxn)
%    maxOrder - maximum order of the TaylorSeries, has to be > abs(intMat) +2
%    potentiation - #nodef
%
% Outputs:
%    scaling - the exponentiation with the chosen algorithm
%
% Example: 
%
% Other m-files required: hornerTaylorSeries.m,taylorSeries.m, intervalMatrixRemainder.m 
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Ivan Brkan
% Written:      23-April-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

potential = 2^potentiation;
if((maxOrder+2)*potential<=norm(intMat,inf))
    scaling= [];
else
    scaling = mpower(hornerTaylorSeries(mtimes(intMat,potential^-1),maxOrder), potential);
end

