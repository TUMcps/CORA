function [p,pTotal]=roadCrossing(simOptions,markovChainSpec)
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

% Author: Matthias Althoff
% Written: 18-June-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%load data
qFree=simOptions.freeDrivingProb;
Psi=simOptions.interactionMatrix;
T=simOptions.transitionMatrix;
pFront=simOptions.frontProbVector;
xSegment=simOptions.xSegment;
field=simOptions.stateField;

%define zero vector z 
z=0*simOptions.initialProbability;

%set up probability distribution of virtual vehicle
pVirtual=z;
pVirtual(xSegment+1)=1;

%get nr of modes
nrOfModes=markovChainSpec.nrOfInputs;

%generate transition filter-------------
%obtain nr of segments
nrOfSeg=get(field,'nrOfSegments');
posFilter=(xSegment-4):xSegment;
tranFilter(prod(nrOfSeg)+1,1)=0;
for iVelocity=1:nrOfSeg(2)
    tranFilter(posFilter+1+(iVelocity-1)*nrOfSeg(1),1)=1;
end
simOptions.tranFilter=tranFilter;
%---------------------------------------

%time step loop
for iStep=1:simOptions.runs
    %initialize mTotal
    mTotal=z;
    
    %input mode loop
    for iMode=1:nrOfModes
        
        %compute m vector for virtual car that is in mode 1
        mVirt(:,iMode)=Psi{iMode,1}*pVirtual;
        
        %if there is another vehicle ahead
        if ~isempty(pFront)
            %initialize m vector
            m(:,iMode)=z;
            %mode loop for leading vehicle
            for iAheadMode=1:nrOfModes
                m(:,iMode)=m(:,iMode)+Psi{iMode,iAheadMode}*pFront.T{iStep}(:,iAheadMode); 
            end  
            %obtain cross selection vector 
            mCross{iStep}(:,iMode)=m(:,iMode);

        %there is no other vehicle ahead
        else
            mAppr{iStep}(:,iMode)=mVirt(:,iMode);
            mCross{iStep}(:,iMode)=ones(prod(nrOfSeg)+1,1)*qFree(iMode);
        end
    end
    
    %combine m with mVirt to mAppr
    if ~isempty(pFront)
        mTemp=(m.*mVirt)';
        mAppr{iStep}=normalizeMatrix(mTemp)';
    end
end

%set transition data
simOptions.pTrans=[];
simOptions.selectionVector=mAppr;

%execute Markov-Chains
[pAppr,pApprTotal,pTrans]=driving(simOptions,markovChainSpec);  

%change input distribution of pTrans in crossing mode
selMode=floor(nrOfModes/2)+1;
for iStep=1:simOptions.runs
    pTransTotal=sum(pTrans{iStep},2);
    for iMode=1:markovChainSpec.nrOfInputs
        if iMode==selMode
            pTrans{iStep}(:,iMode)=pTransTotal;
        else
            pTrans{iStep}(:,iMode)=0*pTransTotal;
        end
    end
end
        
%change options   
simOptions.initialProbability=z;
simOptions.selectionVector=mCross;
simOptions.pTrans=pTrans;
simOptions.tranProb=[];
        
%execute Markov-Chains
[pCross,pCrossTotal]=driving(simOptions,markovChainSpec);   
    
%combine probabilities
for iStep=1:simOptions.runs
    %total probabilities
    pTotal.T{iStep}=pApprTotal.T{iStep}+pCrossTotal.T{iStep};
    pTotal.OT{iStep}=pApprTotal.OT{iStep}+pCrossTotal.OT{iStep};
    
    %input specific probabilities
    p.T{iStep}=pAppr.T{iStep}+pCross.T{iStep};
    p.OT{iStep}=pAppr.OT{iStep}+pCross.OT{iStep};
end
  
%------------- END OF CODE --------------