function [eZ,eI,zPow,iPow,E] = expmMixed(matZ,r,intermediateOrder,maxOrder)
% expmMixed - operator for the exponential matrix of a matrix zonotope,
%    evaluated dependently. Higher order terms are computed via interval
%    arithmetic.
%
% Syntax:
%    [eZ,eI,zPow,iPow,E] = expmMixed(matZ,r,intermediateOrder,maxOrder)
%
% Inputs:
%    matZ - matZonotope object
%    r - time step size
%    intermediate Order - Taylor series order until computation is 
%                           performed with matrix zonotopes
%    maxOrder - maximum Taylor series order until remainder is computed
%
% Outputs:
%    eZ - matrix zonotope exponential part
%    eI - interval matrix exponential part
%    zPow - ?
%    iPow - ?
%    E - ?
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       13-September-2010 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{matZ,'att','matZonotope'}, ...
                {r,'att','numeric',{'scalar','nonnegative'}}, ...
                {intermediateOrder,'att','numeric',{'scalar','integer','>=',2}}, ...
                {maxOrder,'att','numeric',{'scalar','integer'}}});
    
%compute exact terms
[sq,H] = dependentTerms(matZ*(1/r),r);

%init eZ
eZ = H;

%compute powers
zPow=powers(matZ,intermediateOrder,2,sq);

%add first power for input computations
zPow{1}=matZ;

%compute finite Taylor sum
for i=3:intermediateOrder
    eZ = eZ + zPow{i}*(1/factorial(i));
end

%compute interval part
intMat = intervalMatrix(matZ);
[eI,iPow,E] = expm(intMat, r, maxOrder, intermediateOrder+1, ...
    intMat*intervalMatrix(zPow{intermediateOrder}));

% ------------------------------ END OF CODE ------------------------------
