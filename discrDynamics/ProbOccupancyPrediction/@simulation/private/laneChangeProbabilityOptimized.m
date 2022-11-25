function [pChange,plotInfo]=laneChangeProbabilityOptimized...
    (simOptions,selectionVector,pLeft)
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

% Author:       Matthias Althoff
% Written:      15-October-2009
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%load data
ThetaC=simOptions.interactionMatrix;
pBehind=simOptions.frontProbVector{3};

%set weighting vector for motivation computation
wFrontProb=[0 1 2 3 4 4];
wBehindProb=[0 0 0 1 1 1];

%set values for lane change motivation
b=pi/2*tan(0.07);
epsilon=0;

%set iStep
iStep=length(pLeft.T);

%normalize pLeft
pLeft.T{end}=pLeft.T{end}/sum(pLeft.T{end});

%obtain the time step
iStep=length(pLeft.T);

%compute motivation value for staying on the left lane
motLeft=motivationOptimized(wFrontProb,selectionVector.left{iStep},pLeft.T{iStep},simOptions);

%compute motivation value for changing to the right lane
[motRight]=motivationOptimized(wFrontProb,selectionVector.right{iStep},pLeft.T{iStep},simOptions);

%new following vehicle
pLeftAct.T{1}=pLeft.T{iStep};
selectionVecBehind=ThetaC*pLeft.T{end};
%compute motivation value for vehicle driving behind
motBehind=motivationOptimized(wBehindProb,selectionVecBehind,pBehind.T{iStep},simOptions);
%obtain maximum possible morivation for behind vehicle
motBehindMax=sum(wBehindProb);

%compute lane change probability
pChange=2/pi*atan(b*motRight/motLeft*motBehind/motBehindMax)+epsilon;

%gather data for plots
plotInfo.left=motLeft/sum(wFrontProb);
plotInfo.right=motRight/sum(wFrontProb);
plotInfo.behind=motBehind/sum(wBehindProb);
plotInfo.prob=pChange;


%------------- END OF CODE --------------