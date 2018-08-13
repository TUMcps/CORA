function [mat]=sparseDiag(diagVec)
% sparseDiag - creates a sparse diagonal matrix
%
% Syntax:  
%    [mat]=sparseDiag(diagVec)
%
% Inputs:
%    diagVec - vector of diagonal elements
%
% Outputs:
%    mat -sparse diagonal matrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      01-July-2009
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

nrOfElements=length(diagVec);

%create row and column indices
rowInd=1:nrOfElements;
colInd=1:nrOfElements;

%create sparse matrix
mat=sparse(rowInd,colInd,diagVec,nrOfElements,nrOfElements);


%------------- END OF CODE --------------