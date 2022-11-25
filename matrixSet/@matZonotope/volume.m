function vol = volume(matZ)
% volume - computes the volume of a matrix zonotope by computing the volume
% of the corresponding zonotope
%
% Syntax:  
%    vol = volume(matZ)
%
% Inputs:
%    matZ - matrix zonotope
%
% Outputs:
%    vol - volume
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      24-June-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%conversion to a zonotope
Z = zonotope(matZ);

%compute volume of the zonotope
vol = volume(Z);

%------------- END OF CODE --------------