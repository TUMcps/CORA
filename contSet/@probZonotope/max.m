function res = max(probZ,m)
% max - Computes an overapproximation of the maximum on the m-sigma bound
%       according to Eq. (3) in [1]
%
% Syntax:
%    res = max(probZ,m)
%
% Inputs:
%    probZ - probabilistic zonotope object
%    m - m of the m-sigma bound
%
% Outputs:
%    res - overapproximated maximum value
%
% Example:
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2,2);
%    m = 2;
%    res = max(probZ,2);
%
% References:
%    [1] M. Althoff et al. "Safety assessment for stochastic linear systems 
%        using enclosing hulls of probability density functions", ECC 2009
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: probZonotope

% Authors:       Matthias Althoff
% Written:       22-August-2007
% Last update:   08-September-2009
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%obtain covariance matrix
Sigma=sigma(probZ);

%get dimension
d=dim(probZ);

%compute maximum value
res=1/((2*pi)^(d/2)*det(Sigma)^(1/2))*exp(-0.5*m^2);

% ------------------------------ END OF CODE ------------------------------
