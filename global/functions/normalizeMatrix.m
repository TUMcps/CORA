function [M]=normalizeMatrix(M)
% normalize - normalizes a matrix M, such that its columns sum up to one.
%
% Syntax:  
%    [M] = normalizeMatrix(M)
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

% Author:        Matthias Althoff, Mark Wetzlinger
% Written:       27-June-2008
% Last update:   01-May-2020 (MW, vectorization instead of loop)
% Last revision: ---

%------------- BEGIN CODE --------------

%compute column sum
colSum=sum(M);

%get nonzero rows and columns
ind=find(colSum);

%normalize matrix
M(:,ind) = M(:,ind) ./ colSum(ind);


end

%------------- END OF CODE --------------