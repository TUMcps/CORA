function n = dim(I)
% dim - returns the dimension of the ambient space of an STL interval
%
% Syntax:
%    n = dim(I)
%
% Inputs:
%    I - stlInterval object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example:
%    I = stlInterval(0,1);
%    n = dim(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       09-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% STL intervals are time intervals and thus always have dimension 1
n = 1;

% ------------------------------ END OF CODE ------------------------------
