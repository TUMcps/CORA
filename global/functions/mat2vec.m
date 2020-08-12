function vec = mat2vec(mat)
% mat2vec - Stores entries of a matrix into a vector
%
% Syntax:  
%    vec = mat2vec(mat)
%
% Inputs:
%    mat - numerical matrix
%
% Outputs:
%    vec - numerical vector
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      18-June-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%get columns
cols=length(mat(1,:));

%concartenate columns
vec=[];
for i=1:cols
    vec(end+1:end+length(mat(:,i)),1)=mat(:,i);
end

%------------- END OF CODE --------------
