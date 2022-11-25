function [eZ,eI,zPow,iPow,E] = expmIndMixed(matZ,intermediateOrder,maxOrder)
% expmIndMixed - operator for the exponential matrix of a 
% matrix zonotope, evaluated independently. Higher order terms are computed
% via interval arithmetic.
%
% Syntax:  
%    [eZ,eI,zPow,iPow,E] = expmIndMixed(matZ,intermediateOrder,maxOrder)
%
% Inputs:
%    matZ - matrix zonotope
%    intermediate Order - Taylor series order until computation is 
%                           performed with matrix zonotopes
%    maxOrder - maximum Taylor series order until remainder is computed
%
% Outputs:
%    eZ - matrix zonotope exponential part
%    eI - interval matrix exponential part
%    zPow - #nodef
%    iPow - cell array storing the powers of the matrix:
%           A,A^2,...,A^(intermediateOrder)
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
% Written:      18-June-2010 
% Last update:  05-August-2010
% Last revision:---

%------------- BEGIN CODE --------------

%compute powers
zPow=powers(matZ,intermediateOrder);

%compute finite Taylor series
%initialize matrix zonotope
eZ=matZ^0;

%compute finite Taylor sum
for i=1:intermediateOrder
    eZ = eZ + zPow{i}*(1/factorial(i));
end

%compute interval part
intMat = intervalMatrix(matZ);
[eI,iPow,E] = expmInd(intMat, maxOrder, intermediateOrder+1, intMat*intervalMatrix(zPow{intermediateOrder}));



%------------- END OF CODE --------------