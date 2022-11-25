function [p,pTotal]=vehicleFollowing(simOptions,markovChainSpec)
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

% Author: Matthias Althoff
% Written: 18-June-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%load data
Omega=simOptions.interactionMatrix;
T=simOptions.transitionMatrix;
pFront=simOptions.frontProbVector;

%define zero vector z 
z=0*simOptions.initialProbability;

%create selection vectors
for iStep=1:simOptions.runs
    %input mode loop
    for iMode=1:markovChainSpec.nrOfInputs
        %initialize m vector
        m{iStep}(:,iMode)=z;
        %mode loop for leading vehicle
        for iAheadMode=1:markovChainSpec.nrOfInputs
            m{iStep}(:,iMode)=m{iStep}(:,iMode)...
                +Omega{iMode,iAheadMode}*pFront.T{iStep}(:,iAheadMode); 
        end   
    end   
end
simOptions.selectionVector=m;

%set transition data
simOptions.pTrans=[];
simOptions.tranProb=[];

%execute Markov-Chains
[p,pTotal]=driving(simOptions,markovChainSpec);  


%------------- END OF CODE --------------