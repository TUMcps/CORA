function [eI,eI2,iPow,iPow2,E] = expmMixed(matI,r,intermediateOrder,maxOrder)
% expmIndMixed - dummy function for interval matrices.
%
% Syntax:  
%    eZ = expmInd(matZ,maxOrder)
%
% Inputs:
%    matZ - matrix zonotope
%    maxOrder - Taylor series order until computation with interval arith.
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
    [sq,H] = dependentTerms(matI*(1/r),r);

    %init eI
    eI = H;

    %compute powers
    iPow=powers(matI,intermediateOrder,2,sq);
    
    %add first power for input computatins
    iPow{1}=matI;

    %compute finite Taylor sum
    for i=3:intermediateOrder
        eI = eI + iPow{i}*(1/factorial(i));
    end

    %compute interval part
    [eI2,iPow2,E] = expm(matI, r, maxOrder, intermediateOrder+1, matI*iPow{intermediateOrder});

else
    disp('intermediate order too low');
end


%------------- END OF CODE --------------