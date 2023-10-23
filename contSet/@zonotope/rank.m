function r = rank(Z)
% rank - computes the dimension of the affine hull of a zonotope
%
% Syntax:
%    r = rank(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    r - dimension of the affine hull
%
% Example: 
%    Z = zonotope([1 1 0; 0 0 1]);
%    r = rank(Z)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       06-May-2009
% Last update:   15-September-2019 (rename dim -> rank)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

r = rank(Z.G);

% ------------------------------ END OF CODE ------------------------------
