function M = normalizeMatrix(M)
% normalizeMatrix - normalizes a matrix M, such that its columns sum up to
%    one.
%
% Syntax:
%    M = normalizeMatrix(M)
%
% Inputs:
%    M - matrix
%
% Outputs:
%    M - matrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       27-June-2008
% Last update:   01-May-2020 (MW, vectorization instead of loop)
%                20-April-2023 (TL, logical indexing)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%compute column sum
colSum=sum(M);

%get nonzero rows and columns
ind=colSum ~= 0;

%normalize matrix
M(:,ind) = M(:,ind) ./ colSum(ind);

% ------------------------------ END OF CODE ------------------------------
