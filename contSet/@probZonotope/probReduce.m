function [pZ]=probReduce(pZ)
% probReduce - Reduces the number of single Gaussian distributions to n, where
% n is the dimension of the system
%
% Syntax:  
%    [pZ]=probReduce(pZ)
%
% Inputs:
%    pZ - probabilistic zonotope object
%
% Outputs:
%    pZ - probabilistic zonotope object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 27-August-2007
% Last update: 26-February-2008
% Last revision: ---

%------------- BEGIN CODE --------------

if pZ.gauss~=1
    %get new sigma matrix
    G=pZ.g;
    newSigma=G*G';
else
    newSigma=pZ.cov;
end

%ensure symmetry for numerical stability
newSigma=0.5*(newSigma+newSigma');

%get eigenvalue, eigenvectors of newSigma
[V,W]=eig(newSigma);


%compute new generators
pZ.g=V*sqrt(W);

%------------- END OF CODE --------------