function [weightingMat]=steadyStateWeightingMatrix(tranMat,ssVec)
% steadyStateWeightingMatrix - computes the weighting matrix for a transition
% matrix of a Markov-chain, such that the desired steady state vector is
% reached
%
% Syntax:  
%    [weightingMat]=steadyStateWeightingMatrix(tranMat,ssVec)
%
% Inputs:
%    tranMat - transition matrix
%    ssVec - steady state vector
%
% Outputs:
%    weightingMat - weighting matrix
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

tmp=tranMat*ssVec;
tmp2=1./tmp';

weightingMat=ssVec*tmp2;

%------------- END OF CODE --------------