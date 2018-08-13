function [G]=generators(Sigma)
% generators - Returns the generator matrix of a probabilistic zonotope if
% the covariance matrix Sigma is given
%
% Syntax:  
%    [G]=generators(Sigma)
%
% Inputs:
%    Sigma - covariance matrix
%
% Outputs:
%    G - generator vector
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      26-February-2008
% Last update:  09-September-2009
% Last revision:---

%------------- BEGIN CODE --------------

%ensure symmetry for numerical stability
Sigma=0.5*(Sigma+Sigma');

%get eigenvalue, eigenvectors of Sigma
[V,W]=eig(Sigma);

%compute new generators
G=V*sqrt(W);

%------------- END OF CODE --------------