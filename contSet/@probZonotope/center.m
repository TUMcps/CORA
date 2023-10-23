function c = center(probZ)
% center - Returns the center of a probabilistic zonotope
%
% Syntax:
%    c = center(probZ)
%
% Inputs:
%    probZ - probabilistic zonotope object
%
% Outputs:
%    c - center of the probabilistic zonotope
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2);
%    c = center(probZ)
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

c=probZ.Z(:,1);

% ------------------------------ END OF CODE ------------------------------
