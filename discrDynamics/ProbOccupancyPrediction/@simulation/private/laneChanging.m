function [p,pTotal,lcEvolProb]=laneChanging(simOptions,markovChainSpec)
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
%    lcEvolProb - shows the probability of the lane cahnge phases
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 29-October-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%load data
Omega=simOptions.interactionMatrix;
pFrontLeft=simOptions.frontProbVector{1};
pFrontRight=simOptions.frontProbVector{2};
pBehind=simOptions.frontProbVector{3};
initInputProb=simOptions.initialInputProb;
initialProbability=simOptions.initialProbability;
reachIndices=simOptions.reachIndices;

nrOfInputs=markovChainSpec.nrOfInputs;
runs=simOptions.runs;
lcEvolProb=[];

%define zero vector z 
z=0*initialProbability;

selectionVec.left=selectionVector(simOptions,nrOfInputs,z,Omega,pFrontLeft);
selectionVec.right=selectionVector(simOptions,nrOfInputs,z,Omega,pFrontRight);

%set transition data
simOptions.pTrans=[];
simOptions.tranProb=[];

%allow only a single run
simOptions.runs=1;

%init left and right lane probabilities
%set first partial probabilities
for iMode=1:markovChainSpec.nrOfInputs
    pLeft.T{1}(:,iMode)=initInputProb(iMode).*initialProbability; 
    pRight.T{1}(:,iMode)=z;
end 

for iStep=1:runs
   %compute lane changing probability
   [pChange,plotInfo{iStep}]=laneChangeProbability(simOptions,markovChainSpec,selectionVec,pLeft);
   
   %set selection vector
   simOptions.selectionVector=selectionVec.left;   
   %substract changing probability
   pLeft.T{iStep}=(1-pChange)*pLeft.T{iStep};   
   %compute pLeft
   [pLeftNew,pTotalLeftNew]=singleStepDriving(simOptions,markovChainSpec,pLeft);
   pLeft.T{iStep+1}=pLeftNew.T;
   pLeft.OT{iStep}=pTotalLeftNew.OT;
   
   inputDistrLeft=sum(pLeft.T{iStep+1})
   
   %set selection vector
   simOptions.selectionVector=selectionVec.right; 
   %change pLeft and give free probability distribution
   pTmp=sum(pLeft.T{iStep},2);
   pTmp2=pTmp*sparse(simOptions.freeDrivingProb);
   %add changing probability
   pRight.T{iStep}=pRight.T{iStep}+pChange*pTmp2;   
   %compute pRight
   [pRightNew,pTotalRightNew]=singleStepDriving(simOptions,markovChainSpec,pRight);
   pRight.T{iStep+1}=pRightNew.T;
   pRight.OT{iStep}=pTotalRightNew.OT;
   
   inputDistrRight=sum(pRight.T{iStep+1})
   
   
   totalLaneChangeProb=sum(sum(pRight.T{iStep}))
   
   %compute ratio between probability of keeping the lane or changing the
   %lane
%    ratio.left(iStep)=pChange;
%    ratio.right(iStep)=sum(sum(pChange*pLeft.T{iStep}))/sum(sum(pRight.T{iStep}));
   
   lcEvolProb=laneChangeEvolution(pChange,lcEvolProb,markovChainSpec);
end

p=[];
pTotal.left.OT=pLeft.OT;
pTotal.right.OT=pRight.OT;

figure
motivationPlot(plotInfo);

%lcEvolProb
figure
pTmp=lcEvolProb(7,:);
tField=partition([0,5],10);
plotHisto(tField,pTmp);


%--------------------------------------------------------------------------
function [m]=selectionVector(simOptions,nrOfInputs,z,Omega,pFront)


%create selection vectors
for iStep=1:simOptions.runs    
    %input mode loop
    for iMode=1:nrOfInputs
        %initialize m vector
        m{iStep}(:,iMode)=z;
        %mode loop for leading vehicle
        for iAheadMode=1:nrOfInputs
            m{iStep}(:,iMode)=m{iStep}(:,iMode)...
                +Omega{iMode,iAheadMode}*pFront.T{iStep}(:,iAheadMode); 
        end  
    end   
end



%------------- END OF CODE --------------