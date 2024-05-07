function M = randPoint(matZ,varargin)
% randPoint - samples a random matrix from a matrix zonotope
%
% Syntax:
%    M = randPoint(matZ)
%    M = randPoint(matZ,N)
%    M = randPoint(matZ,N,type)
%
% Inputs:
%    matZ - matZonotope object
%    N - number of samples
%    type - type of sampling ('standard' or 'extreme')
%
% Outputs:
%    M - cell-array of sampled matrices
%
% Example: 
%    C = [0 0; 0 0];
%    G{1} = [1 3; -1 2]; G{2} = [2 0; 1 -1];
%    matZ = matZonotope(C,G);
%
%    M = randPoint(matZ,5);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger, Matthias Althoff, Tobias Ladner
% Written:       03-April-2023
% Last update:   15-April-2024 (TL, faster implementation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values
[N,type] = setDefaultValues({1,'standard'},varargin);

% parse input arguments
inputArgsCheck({{matZ,'att','matZonotope'}, ...
            {N,'att','numeric',{'scalar','positive','integer'}},...
            {type,'str',{'standard','extreme'}}});

% init center of points
M = repmat(matZ.C,1,1,N);

% get betas
[n,m,h] = size(matZ.G);
if strcmp(type,'extreme')
    betas = sign(2*rand(n,m,h,N)-1);
elseif strcmp(type, 'standard')
    betas = 2*rand(n,m,h,N)-1;
end

% compute betas*G
if isempty(matZ.G)
    M = M + reshape(sum(pagemtimes(matZ.G,betas),3),n,m,N);
end

% add center


end

% ------------------------------ END OF CODE ------------------------------
