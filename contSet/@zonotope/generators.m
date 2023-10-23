function G = generators(Z)
% generators - returns the generator matrix of a zonotope
%
% Syntax:
%    G = generators(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    G - generator matrix
%
% Example: 
%    Z = zonotope([1 1 0; 0 0 1]);
%    G = generators(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       28-November-2016 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

G = Z.G;

% ------------------------------ END OF CODE ------------------------------
