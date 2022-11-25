function [p,pTotal]=vehicleFollowingOptimized(simOptions,markovChainSpec)
% vehicleFollowing - simulates a vehicle that follows another vehicle
%
% Syntax:  
%    [p,pTotal]=vehicleFollowing(simOptions,markovChainSpec)
%
% Inputs:
%    simOptions - simulation options
%    markovChainSpec - Markov-Chain specifications
%
% Outputs:
%    p - probability distributions for different inputs, times
%    pTotal - total probabilities of different times
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      14-October-2009
% Last update: 	01-August-2016 (using pFront.T{iStep} instead of pFront.T{iStep+1})
% Last revision:---

%------------- BEGIN CODE --------------

%load data
ThetaC=simOptions.interactionMatrix;
pFront=simOptions.frontProbVector;

%define zero vector z 
z=0*simOptions.initialProbability;

%create selection vectors
for iStep=1:(simOptions.runs)
    m{iStep}=ThetaC*pFront.T{iStep}; 
end
simOptions.selectionVector=m;

%set transition data
simOptions.pTrans=[];
simOptions.tranProb=[];

%execute Markov-Chains
[p,pTotal]=drivingOptimized(simOptions,markovChainSpec);  


%------------- END OF CODE --------------