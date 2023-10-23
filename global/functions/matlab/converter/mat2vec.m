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
% See also: vec2mat

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       18-June-2010 
% Last update:   20-April-2023 (TL, simplified)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% reshape
vec = reshape(mat, [], 1);

% ------------------------------ END OF CODE ------------------------------
