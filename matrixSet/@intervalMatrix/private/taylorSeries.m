function taylor = taylorSeries(intMat, maxOrder)
% returns the approximation of e^intMat using the truncated Taylor series
%    with maxOrder iterations.
%
% Syntax:  
%    val = taylorSeries(intMat,maxOrder)
%
% Inputs:
%    intMat - interval matrix (nxn)
%    maxOrder - maximum order of the TaylorSeries, has to be > abs(intMat)+2
%
% Outputs:
%    taylor - the exponentiation with the truncated Taylor series 
%
% Example: 
%
% Other m-files required: intervalMatrixRemainder.m 
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Ivan Brkan
% Written:      06-April-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%get norm and size of the matrix
% needs to be symbolic otherwise it's impossible to calculate the result,
% since the factorials are to big and the some parts of the series are to
% small - inf and 0 are the result

% alpha = norm(intMat);
% if(isequal(intMat.int.sup,intMat.int.inf))
%     alpha = norm(intMat.int.sup,inf);
% else
    alpha = norm(intMat,inf);
% end
n = length(intMat.int.sup);
% if k is not large enough, the error will be to large

if(maxOrder+2 <= alpha)
    taylor = [];
    return;
end

% preparation for the calculation
taylor = zeros(n);
pow = eye(n);

% the calculation of the Taylor series
for i=1:maxOrder
    taylor = plus(taylor,pow);
    pow = mtimes(i^-1,mtimes(pow,intMat));
end
taylor = plus(taylor,pow);

% Remainder needs to be added on the series to mimizize the error
% watch out, if maxOrder is too low, the result won't be helpful
taylor = plus(taylor, intervalMatrixRemainder(intMat,alpha,maxOrder));



