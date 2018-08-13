function [weightingVec]=steadyStateWeighting(tranMat,ssVec)
% steadyStateWeighting - computes the weighting vector for a transition
% matrix of a Markov-chain, such that the desired steady state vector is
% reached
%
% Syntax:  
%    [weightingVec]=steadyStateWeighting(tranMat,ssVec)
%
% Inputs:
%    tranMat - transition matrix
%    ssVec - steady state vector
%
% Outputs:
%    weightingVec - weighting vector
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 25-November-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%obtain number of dimensions
dim=length(ssVec);

for iDim=1:dim
    weightingVec(iDim,1)=ssVec(iDim)/(tranMat(iDim,:)*ssVec);
end

%------------- END OF CODE --------------