function horner = hornerTaylorSeries(intMat,maxOrder) 
% hornerTaylorSeries - returns the approximation of e^intMat using the
%    horner scheme evaluation of the taylorSeries and the 
%
% Syntax:  
%    horner = hornerTaylorSeries(intMat,maxOrder)
%
% Inputs:
%    intMat - intervalMatrix object (nxn)
%    maxOrder - maximum order of the TaylorSeries, has to be greater than
%               abs(intMat) + 2
%
% Outputs:
%    horner - the exponentiation with hornerScheme evaluation of the
%             truncated Taylor series 
%
% Example: 
%    -
%
% Other m-files required: intervalMatrixRemainder.m 
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Ivan Brkan
% Written:      14-April-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%get norm and size of the matrix
alpha = norm(intMat,Inf);

n = length(intMat.int.sup);
% maxOrder has to be alpha +2 in order to calculate a solution
if(maxOrder + 2 <= alpha)
    horner = [];
    return;
end

E = eye(n);
% the algorithm calculates the solution starting at the innermost term,
% so it uses maxOrder first and finishes with 1, as factor the startpoint 
horner = plus(E,mtimes(intMat,maxOrder^-1));
for i= 1:maxOrder-1
    horner = plus(E, mtimes((maxOrder-i)^-1, mtimes(intMat,horner)) );
end

% Remainder has to be added on the solution, in order to minimize the error
horner = plus(horner, intervalMatrixRemainder(intMat,alpha,maxOrder));

%------------- END OF CODE --------------