function [eP,eI] = expmIndMixed(matP,intermediateOrder,maxOrder)
% expmIndMixed - operator for the exponential matrix of a matrix polytope,
%    evaluated independently; higher-order terms are computed via interval
%    arithmetic.
%
% Syntax:
%    [eP,eI] = expmIndMixed(matP,intermediateOrder,maxOrder)
%
% Inputs:
%    matP - matPolytope object
%    intermediateOrder - Taylor series order until computation with interval arith.
%    maxOrder - maximum Taylor series order until remainder is computed
%
% Outputs:
%    eP - matrix polytope exponential part
%    eI - interval matrix exponential part
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       02-July-2010 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%compute finite Taylor series
%initialize matrix polytope
eP=matP^0;
%initialize power
ePpow=matP^0;

%compute finite Taylor sum
for i=1:intermediateOrder
    ePpow = ePpow*matP;
    eP = simplePlus(eP, ePpow*(1/factorial(i)));
end

%compute interval part
intMat = intervalMatrix(matP);
eI = expmInd(intMat, maxOrder, intermediateOrder+1, intMat*intervalMatrix(ePpow));

% ------------------------------ END OF CODE ------------------------------
