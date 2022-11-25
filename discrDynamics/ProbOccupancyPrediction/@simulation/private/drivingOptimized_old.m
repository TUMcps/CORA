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

% Author: Matthias Althoff
% Written: 15-June-2009
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%load data
m=simOptions.selectionVector;
T=simOptions.transitionMatrix;
pAdd=simOptions.pTrans;
tranProb=simOptions.tranProb;
initInputProb=simOptions.initialInputProb;
Gamma=simOptions.Gamma;
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

%compute projection matrix for total probabilities
nrOfInputs=markovChainSpec.nrOfInputs;
nrOfStates=length(simOptions.initialProbability);
nrOfCombinedStates=length(p.T{1});
projMat=sparse(matrixbuilder(nrOfInputs,nrOfStates,1));
projMat(:,1)=[];

%obtain Gfull


%time step loop
for iStep=1:simOptions.runs
             
    %time step solution
    p.T{iStep+1}=T.T*p.T{iStep};
    %time interval solution
    p.OT{iStep}=T.OT*p.T{iStep};

    %compute total probabilities 
    pTotal.T{iStep+1}=projMat*p.T{iStep+1};
    pTotal.OT{iStep}=projMat*p.OT{iStep}; 
    
    %combine m with speed restriction
    m{iStep}=min(m{iStep},speedRes);    
    %combine with free driving probability
    freeDriving=ones(length(m{iStep}),1)*freeDrivingProb;
    m{iStep}=combine(m{iStep},freeDriving);
    
    %compute input transitions
    %get nonzero elements of p.T
    ind=find(pTotal.T{iStep+1});
    %init Gfull
    Gfull=sparse(nrOfCombinedStates,nrOfCombinedStates);
    %obtain input transitions for each state space cell
    for iCell=1:length(ind)
        %obtain index
        index=ind(iCell);
        %compute Gamma matrix
        G=diag(m{iStep}(index,:))*Gamma;
        G=normalizeMatrix(G);

        fullIndices=(1:nrOfInputs)+nrOfInputs*(index-1);
        Gfull(fullIndices,fullIndices)=G;
        Gfull(fullIndices,fullIndices)=sparse(G);
    end
    
    %execute Gamma matrix
    p.T{iStep+1}=Gfull*p.T{iStep+1};    

    
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