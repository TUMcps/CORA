function sig = sigma(probZ)
% sigma - returns Sigma matrix of a probabilistic zonotope
%
% Syntax:
%    sig = sigma(probZ)
%
% Inputs:
%    probZ - probabilistic zonotope object
%
% Outputs:
%    sig - sigma matrix
%
% Example:
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2);
%    sig = sigma(probZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       28-August-2007
% Last update:   26-February-2008
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%reduce probabilistic zonotope first
probZ = probReduce(probZ);

%get new sigma matrix
G = probZ.g;
sig = G*G';

% ------------------------------ END OF CODE ------------------------------
