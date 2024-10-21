function E = exponentialRemainder(intMat,maxOrder)
% exponentialRemainder - returns the remainder of the exponential matrix
%
% Syntax:
%    E = exponentialRemainder(intMat,maxOrder)
%
% Inputs:
%    intMat - intervalMatrix object
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
% See also: 

% Authors:       Matthias Althoff
% Written:       18-June-2010 
% Last update:   14-November-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%compute absolute value bound
M = abs(intMat);

%compute exponential matrix
eM = expm(M);

% no value is infinity
if ~any(any(isnan(eM)))

    %compute first Taylor terms
    Mpow = eye(dim(intMat));
    eMpartial = eye(dim(intMat));
    for i=1:maxOrder
        Mpow = M*Mpow;
        eMpartial = eMpartial + Mpow/factorial(i);
    end

    W = eM-eMpartial;

    %instantiate remainder
    E = intervalMatrix(zeros(dim(intMat)),W);
else
    %instantiate remainder
    E = intervalMatrix(zeros(dim(intMat)),Inf(dim(intMat)));
end

% ------------------------------ END OF CODE ------------------------------
