function [p,pTotal,pTrans]=drivingOptimized_noInteraction(simOptions,markovChainSpec)
% drivingOptimized - simulates a vehicle where the selection matrix has been
% specified; this is an optimized version which is vectorized and thus fits 
% better in the matlab framework.
%
% Syntax:  
%    [p,pTotal]=drivingOptimized(simOptions,markovChainSpec)
%
% Inputs:
%    simOptions - simulation options
%    markovChainSpec - Markov-Chain specifications
%
% Outputs:
%    p - probability distributions for different inputs, times
%    pTotal - total probabilities of different times
%    pTrans - transition probabilities of different times
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      15-June-2009
% Last update:  14-October-2009
% Last revision:---

%------------- BEGIN CODE --------------

%load data
mFull=simOptions.selectionVector;
T=simOptions.transitionMatrix;
initInputProb=simOptions.initialInputProb;
projMat=simOptions.projMat;
GammaFull=simOptions.GammaFull;
speedRes=simOptions.speedRes;
freeDrivingProb=simOptions.freeDrivingProb;

%initialize the total probability
pTotal.T{1}=simOptions.initialProbability;

%define zero vector z 
z=0*pTotal.T{1};

%set first partial probabilities
for iMode=1:markovChainSpec.nrOfInputs
    p.T{1}(:,iMode)=initInputProb(iMode).*pTotal.T{1};       
end 

%reshape initial probability distribution
pTmp=reshape(p.T{1}',prod(size(p.T{1})),1);
p.T{1}=[];
p.T{1}=pTmp;

%get number of inputs and states
nrOfInputs=markovChainSpec.nrOfInputs;
nrOfStates=length(simOptions.initialProbability);
nrOfCombinedStates=length(p.T{1});

%reshape free driving probability distribution
freeDriving=ones(nrOfStates,1)*freeDrivingProb;
freeDrivingTmp=reshape(freeDriving',prod(size(freeDriving)),1);
freeDriving=[];
freeDriving=freeDrivingTmp;


%compute full Gamma (if there is no interaction)------------------
%initialize mAct to 0
mAct=sparse(nrOfCombinedStates,1);    
%combine m with speed restriction
mAct=min(mFull{1},speedRes);    
%combine with free driving probability
mAct=combineOptimized(mAct,freeDriving,nrOfInputs); 

%cGfull=normalizeMatrix(Gfull);ompute Gamma matrix
Gfull=sparse(diag(mAct))*GammaFull;
Gfull=normalizeMatrix(Gfull);
%-----------------------------------------------------------------


%time step loop
for iStep=1:simOptions.runs
             
    %time step solution
    p.T{iStep+1}=T.T*p.T{iStep};
    %time interval solution
    p.OT{iStep}=T.OT*p.T{iStep};

    %compute total probabilities 
    pTotal.T{iStep+1}=projMat*p.T{iStep+1};
    pTotal.OT{iStep}=projMat*p.OT{iStep}; 
    
%     %remove small probabilities
%     ind=find(p.T{iStep+1}>10/nrOfCombinedStates);
%     pTmp=sparse(nrOfCombinedStates,1);    
%     pTmp(ind)=p.T{iStep+1}(ind);
%     pTmp=pTmp/sum(pTmp);
%     p.T{iStep+1}=pTmp;
    
    %get nonzero indices of the reachable set
    nonZeroIndTmp=find(pTotal.T{iStep+1});
    if length(nonZeroIndTmp)>0
        for i=1:length(nonZeroIndTmp)
            nonZeroIndTotal(((i-1)*nrOfInputs+1):(i*nrOfInputs))=(nonZeroIndTmp(i)-1)*nrOfInputs+(1:nrOfInputs);
        end 

        %execute full Gamma matrix
        p.T{iStep+1}(nonZeroIndTotal)=Gfull(nonZeroIndTotal,nonZeroIndTotal)*p.T{iStep+1}(nonZeroIndTotal);
    end
    
end    

%------------- END OF CODE --------------