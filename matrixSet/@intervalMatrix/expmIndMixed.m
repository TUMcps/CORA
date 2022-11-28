function [eI,eI2,iPow,iPow2,E] = expmIndMixed(intMat,intermediateOrder,maxOrder)
% expmIndMixed - ?
%
% Syntax:  
%    expmIndMixed(intMat,intermediateOrder,maxOrder)
%
% Inputs:
%    intMat - intervalMatrix object
%    intermediate Order - Taylor series order until computation is 
%                           performed with matrix zonotopes
%    maxOrder - maximum Taylor series order until remainder is computed
%
% Outputs:
%    eZ - matrix zonotope exponential part
%    eI - interval matrix exponential part
%    iPow - cell array storing the powers of the matrix:
%           A,A^2,...,A^(intermediateOrder)
%    iPow2 - cell array storing the powers of the matrix:
%            A^(intermediateOrder+1),...,A^(maxOrder)
%    E - interval matrix for the remainder
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      13-September-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compute powers
iPow=powers(intMat,intermediateOrder);

%compute finite Taylor series
%initialize matrix zonotope
eI=intMat^0;

%compute finite Taylor sum
for i=1:intermediateOrder
    eI = eI + iPow{i}*(1/factorial(i));
end

%compute interval part
[eI2,iPow2,E] = expmInd(intMat, maxOrder, ...
    intermediateOrder+1, intMat*iPow{intermediateOrder});

%------------- END OF CODE --------------