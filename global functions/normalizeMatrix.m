function [M]=normalizeMatrix(M)
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

[rows,cols]=size(M);

%compute column sum
colSum=sum(M)';

% %get nonzero rows and columns
% ind=find(M);
% [rowInd,colInd]=ind2sub([rows,cols],ind);
% 
% %normalize matrix
% M(ind)=M(ind)./colSum(colInd);

%get nonzero rows and columns
ind=find(colSum);

%normalize matrix
for i=1:length(ind)
    M(:,ind(i))=M(:,ind(i))/colSum(ind(i));
end

%------------- END OF CODE --------------