function mat = vec2mat(varargin)
% vec2mat - Stores entries of a vector in a matrix
%
% Syntax:
%    mat = vec2mat(varargin)
%
% Inputs:
%    vec - vector
%    cols - number of columns for the matrix
%
% Outputs:
%    mat - matrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mat2vec

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       18-June-2010 
% Last update:   22-June-2010
%                05-October-2010
%                20-April-2023 (TL, simplified)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

vec = varargin{1};

% get cols
if nargin==2
    cols = varargin{2};
elseif nargin==1
    % assume square
    cols = sqrt(length(vec));
end

% reshape
mat = reshape(vec, [], cols);

% ------------------------------ END OF CODE ------------------------------
