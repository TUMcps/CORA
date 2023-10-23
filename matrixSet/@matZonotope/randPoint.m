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

% Authors:       Mark Wetzlinger, Matthias Althoff
% Written:       03-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values
[N,type] = setDefaultValues({1,'standard'},varargin);

% parse input arguments
inputArgsCheck({{matZ,'att','matZonotope'}, ...
            {N,'att','numeric',{'scalar','positive','integer'}},...
            {type,'str',{'standard','extreme'}}});

% generate sample matrices
M = cell(N,1);
for i=1:N
    % initialize random matrix
    M{i} = matZ.center;
    
    % add generator matrices
    for iGen=1:matZ.gens
        if strcmp(type,'extreme')
            M{i} = M{i} + sign(2*rand(1)-1) .* matZ.generator{iGen};
        elseif strcmp(type,'standard')
            M{i} = M{i} + (2*rand(1)-1) .* matZ.generator{iGen};
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
