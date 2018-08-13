function [result] = oldMax(pZ,m)
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

% Author: Matthias Althoff
% Written: 22-August-2007
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%load generator matrix
G=1e-4*pZ.pZ.g;
%get dimension and number of generators
[dim,nrOfGen]=size(G);
%generate zero mean zonotope from generators
Z=zonotope([0*G(:,1),G]);

for iGen=1:nrOfGen
    %find nonzero sigmas
    indices=find(pZ.pZ.variance(:,iGen));
    %compute auxiliary value
    auxValue=sum(pZ.pZ.omega(indices,iGen)./pZ.pZ.variance(indices,iGen));
    %compute factors 
    value(iGen)=auxValue*2*norm(pZ.pZ.g(:,iGen));
end

a=nrOfGen-dim+1;
result=prod(value)*exp(-0.5*m^2)^a/(sqrt(2*pi))^nrOfGen/volume(Z);

%------------- END OF CODE --------------