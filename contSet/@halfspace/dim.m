function n = dim(hs)
% dim - returns the dimension of the ambient space of a halfspace
%
% Syntax:
%    n = dim(hs)
%
% Inputs:
%    hs - halfspace object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example: 
%    h = halfspace([1;1],3);
%    dim(h)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       27-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = length(hs.c);

% ------------------------------ END OF CODE ------------------------------
