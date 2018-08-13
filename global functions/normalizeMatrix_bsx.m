function [M]=normalizeMatrix_bsx(M)
% normalize - normalizes a matrix M, such that its columns sum up to one.
%
% Syntax:  
%    [M]=normalizeMatrix(M)
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

% Author: Matthias Althoff
% Written: 27-June-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------


%compute column sum
colSum=sum(M);

%normalize matrix
M=bsxfun(@rdivide, M, colSum);
    

%------------- END OF CODE --------------