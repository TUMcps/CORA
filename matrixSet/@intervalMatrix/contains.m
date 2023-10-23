function res = contains(intMat,M,varargin)
% contains - returns if an interval matrix contains a given matrix
%
% Syntax:
%    res = contains(intMat,M)
%    res = contains(intMat,M,type)
%    res = contains(intMat,M,type,tol)
%
% Inputs:
%    intMat - intervalMatrix object
%    M - matrix
%
% Outputs:
%    res - true/false
%
% Example: 
%    intMat = intervalMatrix([2 3; 1 2],[1 0; 1 1]);
%    M1 = [2.5 3; 2 1];
%    M2 = [4 3.5; 1 2];
%    res1 = contains(intMat,M1);
%    res2 = contains(intMat,M2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       03-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default values
[type,tol] = setDefaultValues({'exact',1e-12},varargin);

% check input arguments
inputArgsCheck({{intMat,'att','intervalMatrix'},...
                {M,'att','numeric'},...
                {type,'str',{'exact','approx'}},...
                {tol,'att','numeric',{'scalar','nonnegative'}}});

% dimension check
equalDimCheck(intMat,M);

% since all values are independent, we use interval/contains for each column
n = dim(intMat);
res = true;
for j=1:n(2)
    if ~contains_(intMat.int(:,j),M(:,j),type,tol)
        res = false; break
    end
end

% ------------------------------ END OF CODE ------------------------------
