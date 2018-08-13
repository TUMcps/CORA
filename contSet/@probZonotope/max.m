function [result] = max(pZ,m)
% max - Computes an overaproximation of the maximum on the m-sigma bound
%
% Syntax:  
%    [result] = max(pZ,m)
%
% Inputs:
%    pZ - probabilistic zonotope object
%    m - m of the m-sigma bound
%
% Outputs:
%    result - overapproximated maximum value
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      22-August-2007
% Last update:  08-September-2009
% Last revision:---

%------------- BEGIN CODE --------------

%obtain covariance matrix
Sigma=sigma(pZ);

%get dimension
dim=length(pZ.g(:,1));

%compute maximum value
result=1/((2*pi)^(dim/2)*det(Sigma)^(1/2))*exp(-0.5*m^2);

%------------- END OF CODE --------------