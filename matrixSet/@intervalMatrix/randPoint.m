function M = randPoint(intMat,varargin)
% randPoint - samples a random matrix from an interval matrix
%
% Syntax:
%    M = randPoint(intMat)
%    M = randPoint(intMat,N)
%    M = randPoint(intMat,N,type)
%
% Inputs:
%    intMat - intervalMatrix object
%    N - number of samples
%    type - type of sampling ('standard' or 'extreme')
%
% Outputs:
%    M - numeric, dimension (n x m x N)
%
% Example: 
%    intMat = intervalMatrix([2 3; 1 2],[1 0; 1 1]);
%    M = randPoint(intMat,5);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger, Matthias Althoff, Tobias Ladner
% Written:       03-April-2023
% Last update:   26-April-2024 (TL, changed output to be a 3d matrix)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values
[N,type] = setDefaultValues({1,'standard'},varargin);

% parse input arguments
inputArgsCheck({{intMat,'att','intervalMatrix'}, ...
            {N,'att','numeric',{'scalar','positive','integer'}},...
            {type,'str',{'standard','extreme'}}});

% obtain dimensions, lower bound, radius
n = intMat.dim;
c = center(intMat.int);
lb = infimum(intMat.int);
r = rad(intMat);

% generate sample matrices
if strcmp(type,'extreme')
    M = c + 0.5*sign(2*rand([n,N])-1) .* r;
elseif strcmp(type,'standard')
    M = lb + rand([n,N]) .* r;
end

% ------------------------------ END OF CODE ------------------------------
