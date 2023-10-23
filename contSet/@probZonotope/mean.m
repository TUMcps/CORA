function Z = mean(probZ)
% mean - Returns the uncertain mean of a probabilistic zonotope
%
% Syntax:
%    Z = mean(probZ)
%
% Inputs:
%    probZ - probabilistic zonotope object
%
% Outputs:
%    Z - uncertain mean of the probabilistic zonotope probZ
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2,2);
%    Z = mean(probZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       27-September-2007 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

Z=zonotope(probZ.Z);

% ------------------------------ END OF CODE ------------------------------
