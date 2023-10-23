function matZ = reduce(matZ,option,order,filterLength)
% reduce - Reduces the order of a matrix zonotope via conversion to
%    zonotopes
%
% Syntax:
%    matZred = reduce(matZ,option,order)
%
% Inputs:
%    matZ - matZonotope object
%    option - selects the reduction method
%    order - order of reduced zonotope
%    filterLength - ?
%
% Outputs:
%    matZ - reduced matrix zonotope
%
% Example: 
%
% Other m-files required: none
% Subfunctions: 
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       24-June-2010
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%convert matrix zonotope to zonotope
Z=zonotope(matZ);

%reduce zonotope
Zred = reduce(Z, option, order, filterLength);

%convert back to matrix zonotope
matZ = matZonotope(Zred);

% ------------------------------ END OF CODE ------------------------------
