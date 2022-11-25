function [p,pTotal]=freeDrivingOptimized(simOptions,markovChainSpec)
% freeDriving - simulates a vehicle that drives on an empty lane
%
% Syntax:  
%    [p,pTotal]=freeDriving(simOptions,markovChainSpec)
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

%get number of cells
[nrOfCells]=length(simOptions.initialProbability);
nrOfModes=markovChainSpec.nrOfInputs;

%create selection vectors
for iStep=1:simOptions.runs
    simOptions.selectionVector{iStep}=ones(nrOfCells*nrOfModes,1);
end

%set transition data
simOptions.pTrans=[];
simOptions.tranProb=[];

%execute Markov-Chains
[p,pTotal]=drivingOptimized(simOptions,markovChainSpec);

%------------- END OF CODE --------------