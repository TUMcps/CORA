function E = intervalMatrixRemainder(intMatr,matrixNorm,maxOrder)
% function, which calculates the remainder for the exponantiation of an
%    intervalMatrix depending on the order of the series
%
% Syntax:  
%    E = intervalMatrixRemainder(intMatr,maxOrder)
%
% Inputs:
%    intMatrix - intervalMatrix 
%    maxOrder - maximum Taylor series order
%
% Outputs
%    E = matrix or intervalMatrix with the remainder for the
%         exponentiation of an intervalMatrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author:       Ivan Brkan
% Written:      03-April-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% Compute the norm 
% a = norm(intMatr,inf);
% getting E started, depending if intMatrix.sup == intMatrix.inf
% it will be [1;1] or [-1;1]

a = matrixNorm;
n= length(intMatr.int.sup);

if(isequal(intMatr.int.sup,intMatr.int.inf))
    E=intervalMatrix(ones(n),0);
    % a= norm(intMatr.int.sup,inf);
else
    E= intervalMatrix(zeros(n),ones(n)); 
    % a= norm(intMatr,inf);
end

if(a>=maxOrder+2)
    E = [];
    return;
end
% no value is infinity .from exponentialRemainder.m

%if ( ~any(any(isnan(intMatr.int.inf))) && ~any(any(isnan(intMatr.int.sup))) )
    
   %compute fomula p(a,K)= a^(K+1)/((K+1)!(1-a/(K+2)))
   %factor = a^(maxOrder+1) / ( factorial(maxOrder+1)*(1-(a /(maxOrder+2))) );
   %factor = double( sym(a)^(maxOrder+1) / (factorial(sym(maxOrder +1))*(1-(a/(maxOrder+2))) ));
  
    increasedOrder = floor(maxOrder)+1; 
    counter = a*ones(1,(increasedOrder));
    denominator = 1:(increasedOrder);
    tmp = counter ./ denominator;
    %tmp = tmp(~(tmp==0));
    factor = prod(tmp);
    factor = factor / (1-(a/(maxOrder+2)));

    E= E*factor;
    
%else
    %instantiate remainder - from exponentialRemainder.m
    %E = [];
%end

end

