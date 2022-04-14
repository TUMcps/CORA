function G = generators(Z)
% generators - returns the generators of a zonotope as a matrix whose
%    column vectors are the generators
%
% Syntax:  
%    G = generators(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    G - matrix of generators stored as column vectors
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

% Author:       Matthias Althoff
% Written:      28-November-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

G = Z.Z(:,2:end);

%------------- END OF CODE --------------