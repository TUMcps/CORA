function n = dim(O)
% dim - returns the dimension of the ambient space of an empty set
%
% Syntax:
%    n = dim(O)
%
% Inputs:
%    O - emptySet object
%
% Outputs:
%    n - dimension
%
% Example: 
%    O = emptySet(2);
%    n = dim(O);
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

n = O.dimension;

% ------------------------------ END OF CODE ------------------------------
