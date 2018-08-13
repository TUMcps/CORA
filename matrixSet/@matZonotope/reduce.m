function matZred = reduce(matZ,option,order, filterLength)
% reduce - Reduces the order of a matrix zonotope
% This is done by converting the matrix zonotope to a zonotope
%
% Syntax:  
%    matZred = reduce(matZ,option,order)
%
% Inputs:
%    matZ - matrix zonotope
%    option - selects the reduction method
%    order - order of reduced zonotope
%
% Outputs:
%    matZred - reduced matrix zonotope
%
% Example: 
%
% Other m-files required: none
% Subfunctions: 
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      24-June-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%convert matrix zonotope to zonotope
Z=zonotope(matZ);

%reduce zonotope
Zred = reduce(Z, option, order, filterLength);

%convert back to matrix zonotope
matZred = matZonotope(Zred);

%------------- END OF CODE --------------