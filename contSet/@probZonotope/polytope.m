function P = polytope(probZ)
% polytope - Converts the mean of a probabilistic zonotope
%    to a polytope representation
%
% Syntax:
%    P = polytope(probZ)
%
% Inputs:
%    probZ - probabilistic zonotope object
%
% Outputs:
%    P - polytope object
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2,2);
%    P = polytope(probZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  vertices

% Authors:       Matthias Althoff
% Written:       18-September-2007
% Last update:   17-August-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%reduce probabilistic zonotope
P = polytope(zonotope(probZ.Z));

% ------------------------------ END OF CODE ------------------------------
