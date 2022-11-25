function [p,pTotal]=roadCrossingOptimized(simOptions,markovChainSpec)
% roadCrossing - simulates a vehicle crossing an intersection that has no
% right of way; if simOptions.otherProbVectorSeries=[], no vehicle
% driving ahead is considered, otherwise it is assumed that this
% probability distribution refers to the car driving ahead
%
% Syntax:  
%    [p,pTotal]=roadCrossing(simOptions,markovChainSpec)
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
% Last update:  13-August-2018
% Last revision:---

%------------- BEGIN CODE --------------

%load data
ThetaC=simOptions.interactionMatrix;
pFront=simOptions.frontProbVector;
xSegment=simOptions.xSegment;
field=simOptions.stateField;

%get number of cells
nrOfCells=length(simOptions.initialProbability);
nrOfModes=markovChainSpec.nrOfInputs;

%set up probability distribution of virtual vehicle
pVirtual=zeros(nrOfCells*nrOfModes,1);
pVirtual(xSegment*nrOfModes+1)=1;


%generate transition filter-------------
%obtain nr of segments
nrOfSeg=field.nrOfSegments;
posFilter=(xSegment-6):xSegment;
tranFilter(nrOfCells*nrOfModes,1)=0;
indices=[];
for iVelocity=1:nrOfSeg(2)
    stateIndices=posFilter+1+(iVelocity-1)*nrOfSeg(1);
    indicesTmp=(stateIndices-1)*nrOfModes;
    for i=1:length(indicesTmp)
        indices(end+1:end+nrOfModes)=indicesTmp(i)+(1:nrOfModes);
    end
    tranFilter(indices,1)=1;
end
simOptions.tranFilter=sparse(tranFilter);
%---------------------------------------

%time step loop
%for iStep=1:(simOptions.runs-1)
for iStep=1:(simOptions.runs)
    
    %compute m vector for virtual car
    mVirt=ThetaC*pVirtual;
        
    %if there is another vehicle ahead
    if ~isempty(pFront)
        m=ThetaC*pFront.T{iStep}; 

    %there is no other vehicle ahead
    else
        m=ones(nrOfCells*nrOfModes,1);
    end
    
    %combine m with mVirt to mAppr (approach)
    mAppr{iStep}=min(mVirt,m);
    mCross{iStep}=m;
end

%set transition data
simOptions.pTrans=[];
simOptions.selectionVector=mAppr;

%execute Markov-Chains for approaching phase
[pAppr,pApprTotal,pTrans]=drivingOptimized(simOptions,markovChainSpec);  

%change input distribution of pTrans in crossing mode
selMode=floor(nrOfModes/2)+1;
for iStep=1:(simOptions.runs)
    pTransTotal=simOptions.projMat*pTrans{iStep+1};
    for iMode=1:markovChainSpec.nrOfInputs
        if iMode==selMode
            pTransTmp(:,iMode)=pTransTotal;
        else
            pTransTmp(:,iMode)=0*pTransTotal;
        end
    end
    %bring transition probabilities to new form
    pTmp=reshape(pTransTmp',prod(size(pTransTmp)),1);
    pTrans{iStep+1}=[];
    pTrans{iStep+1}=pTmp;
end
        
%change options for crossing phase
simOptions.initialProbability=zeros(nrOfCells,1);
simOptions.selectionVector=mCross;
simOptions.pTrans=pTrans;
simOptions.tranProb=[];
        
%execute Markov-Chains
[pCross,pCrossTotal]=drivingOptimized(simOptions,markovChainSpec);   
    
%combine probabilities
for iStep=1:(simOptions.runs)
    %total probabilities
    pTotal.T{iStep}=pApprTotal.T{iStep}+pCrossTotal.T{iStep};
    pTotal.OT{iStep}=pApprTotal.OT{iStep}+pCrossTotal.OT{iStep};
    
    %input specific probabilities
    p.T{iStep}=pAppr.T{iStep}+pCross.T{iStep};
    p.OT{iStep}=pAppr.OT{iStep}+pCross.OT{iStep};
end
  
%------------- END OF CODE --------------