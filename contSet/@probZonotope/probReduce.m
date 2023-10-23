function probZ = probReduce(probZ)
% probReduce - Reduces the number of single Gaussian distributions to
%    the dimension of the system
%
% Syntax:
%    probZ = probReduce(probZ)
%
% Inputs:
%    probZ - probabilistic zonotope object
%
% Outputs:
%    probZ - probabilistic zonotope object
%
% Example:
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2);
%    probZred = probReduce(probZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       27-August-2007
% Last update:   26-February-2008
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if probZ.gauss~=1
    %get new sigma matrix
    G=probZ.g;
    newSigma=G*G';
else
    newSigma=probZ.cov;
end

%ensure symmetry for numerical stability
newSigma=0.5*(newSigma+newSigma');

%get eigenvalue, eigenvectors of newSigma
[V,W]=eig(newSigma);

%compute new generators
probZ.g=V*sqrt(W);

% ------------------------------ END OF CODE ------------------------------
