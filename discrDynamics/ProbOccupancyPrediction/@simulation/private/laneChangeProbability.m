function [pChange,plotInfo]=laneChangeProbability...
    (simOptions,markovChainSpec,selectionVector,pLeft)
% laneChanging - simulates a vehicle changing lanes on a multi-lane road
%
% Syntax:  
%    [p,pTotal]=laneChanging(simOptions,markovChainSpec)
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
% Written: 30-October-2008
% Last update: 27-March-2009
% Last revision: ---

%------------- BEGIN CODE --------------

%load data
Omega=simOptions.interactionMatrix;
pBehind=simOptions.frontProbVector{3};
reachIndices=simOptions.reachIndices;

nrOfInputs=markovChainSpec.nrOfInputs;
runs=simOptions.runs;

%set weighting vector for motivation computation
wFrontProb=[0 1 2 3 3];
wBehindProb=[0 0 1 1 1];

%set values for lane change motivation
b=pi/2*tan(0.07);
epsilon=0;

%define zero vector z 
z=0*simOptions.initialProbability;

%set iStep
iStep=length(pLeft.T);

%set transition data
simOptions.pTrans=[];
simOptions.tranProb=[];


%get sum of the probabilities
s=sum(sum(pLeft.T{end}));
pLeft.T{end}=pLeft.T{end}/s;

%obtain the time step
iStep=length(pLeft.T);

%compute motivation value for staying on the left lane
motLeft=motivation(wFrontProb,selectionVector.left{iStep},pLeft.T{iStep},simOptions);

%compute motivation value for changing to the right lane
[motRight]=motivation(wFrontProb,selectionVector.right{iStep},pLeft.T{iStep},simOptions);

%NEW FOLLOWING VEHICLE
pLeftAct.T{1}=pLeft.T{iStep};
selectionVecBehind=singleStepSelectionVector(simOptions,nrOfInputs,z,Omega,pLeft.T{end});
%compute motivation value for vehicle driving behind
motBehind=motivation(wBehindProb,selectionVecBehind,pBehind.T{iStep},simOptions);
%obtain maximum possible morivation for behind vehicle
motBehindMax=sum(wBehindProb);

%COMPUTE LANE CHANGE PROBABILITY
pChange=2/pi*atan(b*motRight/motLeft*motBehind/motBehindMax)+epsilon;


plotInfo.left=motLeft/sum(wFrontProb);
plotInfo.right=motRight/sum(wFrontProb);
plotInfo.behind=motBehind/sum(wBehindProb);
plotInfo.prob=pChange;

%--------------------------------------------------------------------------
function [m]=singleStepSelectionVector(simOptions,nrOfInputs,z,Omega,pFront)


%input mode loop
for iMode=1:nrOfInputs
    %initialize m vector
    m(:,iMode)=z;
    %mode loop for leading vehicle
    for iAheadMode=1:nrOfInputs
        m(:,iMode)=m(:,iMode)...
            +Omega{iMode,iAheadMode}*pFront(:,iAheadMode); 
    end  
end   


%------------- END OF CODE --------------