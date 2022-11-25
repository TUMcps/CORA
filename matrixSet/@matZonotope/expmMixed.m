function [eZ,eI,zPow,iPow,E] = expmMixed(matZ,r,intermediateOrder,maxOrder)
% expmMixed - operator for the exponential matrix of a 
% matrix zonotope, evaluated dependently. Higher order terms are computed
% via interval arithmetic.
%
% Syntax:  
%    [eZ,eI,zPow,iPow,E] = expmMixed(matZ,r,intermediateOrder,maxOrder)
%
% Inputs:
%    matZ - matrix zonotope
%    r - time increment
%    intermediate Order - Taylor series order until computation is 
%                           performed with matrix zonotopes
%    maxOrder - maximum Taylor series order until remainder is computed
%
% Outputs:
%    eZ - matrix zonotope exponential part
%    eI - interval matrix exponential part
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

if (intermediateOrder>=2)
    
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
    matI = intervalMatrix(matZ);
    [eI,iPow,E] = expm(matI, r, maxOrder, intermediateOrder+1, matI*intervalMatrix(zPow{intermediateOrder}));

else
    disp('intermediate order too low');
end


%------------- END OF CODE --------------