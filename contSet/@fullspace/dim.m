function n = dim(fs)
% dim - returns the dimension of the ambient space of a full-dimensional
%    space
%    case R^0: 0
%
% Syntax:
%    n = dim(fs)
%
% Inputs:
%    fs - fullspace object
%
% Outputs:
%    n - dimension
%
% Example: 
%    fs = fullspace(2);
%    n = dim(fs);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = fs.dimension;

% ------------------------------ END OF CODE ------------------------------
