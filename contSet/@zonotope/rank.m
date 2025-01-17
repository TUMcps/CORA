function r = rank(Z,varargin)
% rank - computes the dimension of the affine hull of a zonotope
%
% Syntax:
%    r = rank(Z)
%
% Inputs:
%    Z - zonotope object
%    tol - numeric, tolerance
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
%                15-January-2024 (TL, added tol)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
tol = setDefaultValues({0},varargin);

% compute rank
r = rank(Z.G, tol);

% ------------------------------ END OF CODE ------------------------------
