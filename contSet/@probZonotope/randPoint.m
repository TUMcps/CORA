function p = randPoint(probZ)
% randPoint - generates a random point of a probabilistic 
% zonotope
%
% Syntax:  
%    p = randPoint(probZ)
%
% Inputs:
%    probZ - probabilistic zonotope object
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2);
%    p = randPoint(probZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: randPoint

% Author:       Matthias Althoff
% Written:      17-July-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% obtain extreme point of zonotope
c = randPoint(zonotope(probZ.Z));

% generate random point from normal distribution
p = mvnrnd(c,sigma(probZ))';

%------------- END OF CODE --------------