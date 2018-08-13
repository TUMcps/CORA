function E = exponentialRemainder(intMat,maxOrder)
% exponentialRemainder - returns the remainder of the exponential matrix
%
% Syntax:  
%    E = exponentialRemainder(intMat,maxOrder)
%
% Inputs:
%    intMat - interval matrix
%    maxOrder - maximum order of Taylor series
%
% Outputs:
%    E - remainder of exponential 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author:       Matthias Althoff
% Written:      18-June-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compute absolute value bound
M = abs(intMat);

%compute exponential matrix
eM = expm(M);

%compute first Taylor terms
Mpow = eye(intMat.dim);
eMpartial = eye(intMat.dim);
for i=1:maxOrder
    Mpow = M*Mpow;
    eMpartial = eMpartial + Mpow/factorial(i);
end

W = eM-eMpartial;

%instantiate remainder
E = intervalMatrix(zeros(intMat.dim),W);

%------------- END OF CODE --------------