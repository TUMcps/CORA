function [eI,eI2,iPow,iPow2,E] = expmMixed(intMat,r,intermediateOrder,maxOrder)
% expmMixed - ?
%
% Syntax:  
%    [eI,eI2,iPow,iPow2,E] = expmMixed(intMat,r,intermediateOrder,maxOrder)
%
% Inputs:
%    intMat - interval matrix object
%    r - time step size
%    intermediateOrder - max taylor order for the first part of evaluation
%    maxOrder - max taylor order for the second part of evaluation
%
% Outputs:
%    eI - interval matrix for the first part of evaluation
%    eI2 - interval matrix for the second part of evaluation
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

if intermediateOrder >= 2
    
    %compute exact terms
    [sq,H] = dependentTerms(intMat*(1/r),r);

    %init eI
    eI = H;

    %compute powers
    iPow=powers(intMat,intermediateOrder,2,sq);
    
    %add first power for input computatins
    iPow{1}=intMat;

    %compute finite Taylor sum
    for i=3:intermediateOrder
        eI = eI + iPow{i}*(1/factorial(i));
    end

    %compute interval part
    [eI2,iPow2,E] = expm(intMat, r, ...
        maxOrder, intermediateOrder+1, intMat*iPow{intermediateOrder});

else
    throw(CORAerror('CORA:wrongValue','third','Intermediate order too low'));
end


%------------- END OF CODE --------------