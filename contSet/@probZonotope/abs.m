function probZ = abs(probZ)
% abs - Returns a probabilistic zonotope with absolute values
%    of the center and the interval generators
%
% Syntax:
%    probZ = abs(probZ)
%
% Inputs:
%    probZ - probZonotope object
%
% Outputs:
%    probZ - probZonotope, where Z=(|c|,|g_1|,...,|g_n|)
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2);
%    probZabs = abs(probZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       30-September-2006 
% Last update:   22-March-2007
%                10-June-2020 (MW, update header)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

probZ.Z = abs(probZ.Z);

% ------------------------------ END OF CODE ------------------------------
