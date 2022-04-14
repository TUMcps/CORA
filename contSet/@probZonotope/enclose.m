function probZ = enclose(probZ,Ar,varargin)
% enclose - Generates a probabilistic zonotope that encloses two 
%    probabilistic zonotopes probZ, A*probZ according to Sec. VI.A in [1]
%
% Syntax:  
%    probZ = enclose(pZ,A)
%
% Inputs:
%    probZ - first probabilistic zonotope object
%    Ar - system matrix multiplied with time increment r
%    Rtrans - (optional)
%
% Outputs:
%    probZ - probabilistic zonotope enclosing pZ and A*pZ
%
% Example:
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2);
%    A = rand(2,2);
%    enclose(probZ,A);
%
% References:
%    [1] M. Althoff et al. "Safety assessment for stochastic linear systems 
%        using enclosing hulls of probability density functions", ECC 2009
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      06-September-2007
% Last update:  03-September-2009
%               04-September-2009
%               17-July-2020
% Last revision:---

%------------- BEGIN CODE --------------

if nargin < 3
    Rtrans = zeros(length(Ar),1);
else
    Rtrans = varargin{1};
end

%TO DO: change computation as it is in the dissertation!

%get dimension
d = dim(probZ);

%retrieve uncertain zonotope deltaZ
Zaux = (expm(Ar)-eye(d))*zonotope(probZ)+Rtrans;
deltaZ = enclose(zonotope(zeros(d,1)),Zaux);
probZ = probZ+deltaZ;

%change probabilistic generators if det(A)<1
if trace(Ar)<0
    if probZ.gauss
        probZ.cov=expm(Ar)*probZ.cov*expm(Ar)';
    else
        probZ.g=expm(Ar)*probZ.g;
    end
end

%------------- END OF CODE --------------