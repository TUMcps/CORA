function [p,pTotal,lcEvolProb]=laneChangingOptimized(simOptions,markovChainSpec)
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

% Author:       Matthias Althoff
% Written:      15-October-2009
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%load data
ThetaC=simOptions.interactionMatrix;
projMat=simOptions.projMat;
pFrontLeft=simOptions.frontProbVector{1};
pFrontRight=simOptions.frontProbVector{2};
pBehind=simOptions.frontProbVector{3};
initInputProb=simOptions.initialInputProb;
initialProbability=simOptions.initialProbability;

%nr of inputs, states and runs
nrOfInputs=markovChainSpec.nrOfInputs;
nrOfStates=length(simOptions.initialProbability);
runs=simOptions.runs;

%generate input projection matrix
inputProjMat=matrixbuilder(nrOfInputs,nrOfStates,0);
inputProjMat(:,1)=[];
simOptions.inputProjMat=inputProjMat;

lcEvolProb=[];

selectionVec.left=selectionVector(runs,ThetaC,pFrontLeft);
selectionVec.right=selectionVector(runs,ThetaC,pFrontRight);

%set transition data
simOptions.pTrans=[];
simOptions.tranProb=[];

%allow only a single run
simOptions.runs=1;

%init left and right lane probabilities
for iMode=1:markovChainSpec.nrOfInputs
    pLeft.T{1}(:,iMode)=initInputProb(iMode).*initialProbability; 
    pRight.T{1}(:,iMode)=0*initialProbability;
end 
%reshape initial probability distribution
pLeft.T{1}=re_shape(pLeft.T{1});
pRight.T{1}=re_shape(pRight.T{1});



for iStep=1:runs
   %compute lane changing probability
   [pChange,plotInfo{iStep}]=laneChangeProbabilityOptimized(simOptions,selectionVec,pLeft);
   
   %set selection vector
   simOptions.selectionVector=selectionVec.left;   
   %substract changing probability
   pLeft.T{iStep}=(1-pChange)*pLeft.T{iStep};   
   %compute pLeft
   [pLeftNew,pTotalLeftNew]=singleStepDrivingOptimized(simOptions,pLeft);
   pLeft.T{iStep+1}=pLeftNew.T;
   pLeft.OT{iStep}=pLeftNew.OT;
   pTotalLeft.T{iStep+1}=pTotalLeftNew.T;
   pTotalLeft.OT{iStep}=pTotalLeftNew.OT;
   
   inputDistrLeft=inputProjMat*pLeft.T{iStep+1};
   
   %set selection vector
   simOptions.selectionVector=selectionVec.right; 
   %change pLeft and give free probability distribution
   pTmp=projMat*pLeft.T{iStep};
   pTmp2=pTmp*sparse(simOptions.freeDrivingProb);
   pTmp2=re_shape(pTmp2);
   
   %add changing probability
   pRight.T{iStep}=pRight.T{iStep}+pChange*pTmp2;   
   %compute pRight
   [pRightNew,pTotalRightNew]=singleStepDrivingOptimized(simOptions,pRight);
   pRight.T{iStep+1}=pRightNew.T;
   pRight.OT{iStep}=pRightNew.OT;
   pTotalRight.T{iStep+1}=pTotalRightNew.T;
   pTotalRight.OT{iStep}=pTotalRightNew.OT;
   
   inputDistrRight=inputProjMat*pRight.T{iStep+1};
   
   
   totalLaneChangeProb=sum(pRight.T{iStep})
   
   %compute ratio between probability of keeping the lane or changing the
   %lane
%    ratio.left(iStep)=pChange;
%    ratio.right(iStep)=sum(sum(pChange*pLeft.T{iStep}))/sum(sum(pRight.T{iStep}));
   
   lcEvolProb=laneChangeEvolution(pChange,lcEvolProb,markovChainSpec);
end

p.left.T=pLeft.T;
p.left.OT=pLeft.OT;
p.right.T=pRight.T;
p.right.OT=pRight.OT;

pTotal.left.OT=pTotalLeft.OT;
pTotal.right.OT=pTotalRight.OT;

% figure
% motivationPlot(plotInfo);
% 
% %lcEvolProb
% figure
% pTmp=lcEvolProb(7,:);
% tField=partition([0,5],10);
% plotHisto(tField,pTmp);


%--------------------------------------------------------------------------
function [m]=selectionVector(runs,ThetaC,pFront)
%create selection vectors
for iStep=1:runs  
    m{iStep}=ThetaC*pFront.T{iStep};
end

%--------------------------------------------------------------------------
function [pTmp]=re_shape(p)

pTmp=reshape(p',prod(size(p)),1);

%------------- END OF CODE --------------