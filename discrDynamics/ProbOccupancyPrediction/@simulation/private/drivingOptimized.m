function [p,pTotal,pTrans]=drivingOptimized(simOptions,markovChainSpec)
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
% Last update:  01-August-2016
%               14-August-2018
% Last revision:---

%------------- BEGIN CODE --------------

%load data
mFull=simOptions.selectionVector;
T=simOptions.transitionMatrix;
pAdd=simOptions.pTrans;
tranProb=simOptions.tranProb;
initInputProb=simOptions.initialInputProb;
projMat=simOptions.projMat;
GammaFull=simOptions.GammaFull;
speedRes=simOptions.speedRes;
freeDrivingProb=simOptions.freeDrivingProb;

%initialize the total probability
pTotal.T{1}=simOptions.initialProbability;


%set first partial probabilities
for iMode=1:markovChainSpec.nrOfInputs
    p.T{1}(:,iMode)=initInputProb(iMode).*pTotal.T{1};       
end 

%reshape initial probability distribution
pTmp=reshape(p.T{1}',numel(p.T{1}),1);
p.T{1}=[];
p.T{1}=pTmp;

%get number of inputs and states
nrOfInputs=markovChainSpec.nrOfInputs;
nrOfStates=length(simOptions.initialProbability);
nrOfCombinedStates=length(p.T{1});

%reshape free driving probability distribution
freeDriving=ones(nrOfStates,1)*freeDrivingProb;
freeDrivingTmp=reshape(freeDriving',numel(freeDriving),1);
freeDriving=freeDrivingTmp;


%time step loop
for iStep=1:(simOptions.runs) 
             
    %time step solution
    p.T{iStep+1}=T.T*p.T{iStep};
    %time interval solution
    p.OT{iStep}=T.OT*p.T{iStep};

    %compute total probabilities 
    pTotal.T{iStep+1}=projMat*p.T{iStep+1};
    pTotal.OT{iStep}=projMat*p.OT{iStep}; 
    
    %remove small probabilities
    ind=find(p.T{iStep+1}>10/nrOfCombinedStates);
    %ind=find(p.T{iStep+1}>3/nrOfCombinedStates);
    pTmp=sparse(nrOfCombinedStates,1);    
    pTmp(ind)=p.T{iStep+1}(ind);
    p.T{iStep+1}=pTmp;
    
    %get nonzero indices of the reachable set
    nonZeroIndTmp=find(pTotal.T{iStep+1});
    if ~isempty(nonZeroIndTmp)
        for i=1:length(nonZeroIndTmp)
            nonZeroIndTotal(((i-1)*nrOfInputs+1):(i*nrOfInputs))=(nonZeroIndTmp(i)-1)*nrOfInputs+(1:nrOfInputs);
        end

  
        %combine m with speed restriction
        mPartial=min(mFull{iStep}(nonZeroIndTotal),speedRes(nonZeroIndTotal));   

        %combine with free driving probability
        mPartial=combineOptimized(mPartial,freeDriving(nonZeroIndTotal),nrOfInputs); 

        %compute Gamma matrix
        Gpartial=sparse(diag(mPartial))*GammaFull(nonZeroIndTotal,nonZeroIndTotal);
        Gpartial=normalizeMatrix(Gpartial);
        %execute Gamma matrix
        p.T{iStep+1}(nonZeroIndTotal)=Gpartial*p.T{iStep+1}(nonZeroIndTotal);    
    end
    
    %probability outflow
    if ~isempty(tranProb)
        %input mode loop         
        pTrans{iStep+1}=tranProb(iStep)*diag(simOptions.tranFilter)*p.T{iStep+1};
        p.T{iStep+1}=p.T{iStep+1}-pTrans{iStep+1};
    else
        pTrans=[];
    end   
    
    %probability inflow
    if ~isempty(pAdd)
        p.T{iStep+1}=p.T{iStep+1}+pAdd{iStep+1};
    end
end    

%------------- END OF CODE --------------