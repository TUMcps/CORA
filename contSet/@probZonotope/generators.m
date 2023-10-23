function G = generators(probZ)
% generators - Returns the generator matrix of a probabilistic zonotope
%    using its covariance matrix Sigma
%
% Syntax:
%    G = generators(probZ)
%
% Inputs:
%    probZ - probabilistic zonotope
%
% Outputs:
%    G - generator vector
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2,2);
%    G = generators(probZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       26-February-2008
% Last update:   09-September-2009
%                10-June-2020 (MW, rewrite using probZonotope object)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%ensure symmetry for numerical stability
Sigma=0.5*(probZ.cov+probZ.cov');

%get eigenvalue, eigenvectors of Sigma
[V,W]=eig(Sigma);

%compute new generators
G=V*sqrt(W);

% ------------------------------ END OF CODE ------------------------------
