function probZ = enclose(probZ,Ar,varargin)
% enclose - encloses a probabilistic zonotope and its affine transformation
%    according to Sec. VI.A in [1]
%
% Description:
%    Computes the set
%    { a x1 + (1 - a) * x2 | x1 \in probZ, x2 \in Z2, a \in [0,1] }
%    where Z2 = Ar*probZ + v
%
% Syntax:
%    probZ = enclose(probZ,Ar)
%
% Inputs:
%    probZ - probZonotope object
%    Ar - system matrix multiplied with time increment r
%    Rtrans - (optional) point
%
% Outputs:
%    probZ - probabilistic zonotope enclosing pZ and Ar*pZ
%
% Example:
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2);
%    Ar = rand(2,2);
%    enclose(probZ,Ar);
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

% Authors:       Matthias Althoff
% Written:       06-September-2007
% Last update:   03-September-2009
%                04-September-2009
%                17-July-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
Rtrans = setDefaultValues({zeros(length(Ar),1)},varargin);

% check input arguments
inputArgsCheck({{probZ,'att','probZonotope'};
                {Ar,'att','numeric','nonnan'};
                {Rtrans,'att','numeric','nonnan'}});

%TODO: change computation as it is in the dissertation!

%get dimension
n = dim(probZ);

%retrieve uncertain zonotope deltaZ
Zaux = (expm(Ar)-eye(n))*zonotope(probZ)+Rtrans;
deltaZ = enclose(zonotope(zeros(n,1)),Zaux);
probZ = probZ+deltaZ;

%change probabilistic generators if det(A)<1
if trace(Ar)<0
    if probZ.gauss
        probZ.cov=expm(Ar)*probZ.cov*expm(Ar)';
    else
        probZ.g=expm(Ar)*probZ.g;
    end
end

% ------------------------------ END OF CODE ------------------------------
