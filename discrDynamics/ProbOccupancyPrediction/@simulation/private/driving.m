function [p,pTotal,pTrans]=driving(simOptions,markovChainSpec)
% driving - simulates a vehicle where the selection matrix has been
% specified
%
% Syntax:  
%    [p,pTotal]=driving(simOptions,markovChainSpec)
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
% Written: 03-November-2008
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
reachIndices=simOptions.reachIndices;
freeDrivingProb=simOptions.freeDrivingProb;

%initialize the total probability
pTotal.T{1}=simOptions.initialProbability;

%define zero vector z 
z=0*pTotal.T{1};

%set first partial probabilities
for iMode=1:markovChainSpec.nrOfInputs
    p.T{1}(:,iMode)=initInputProb(iMode).*pTotal.T{1};       
end 

%time step loop
for iStep=1:simOptions.runs
    
    %initialize pTotal
    pTotal.T{iStep+1}=z;
    pTotal.OT{iStep}=z;
    
    if iStep<simOptions.runs %<-- Alexanders Korrektur
    %input mode loop
    for iMode=1:markovChainSpec.nrOfInputs     
        %time step solution
        p.T{iStep+1}(:,iMode)=T.T{iMode}*p.T{iStep}(:,iMode);
        %time interval solution
        p.OT{iStep}(:,iMode)=T.OT{iMode}*p.T{iStep}(:,iMode);
        
        %compute total probabilities 
        pTotal.T{iStep+1}=pTotal.T{iStep+1}+p.T{iStep+1}(:,iMode);
        pTotal.OT{iStep}=pTotal.OT{iStep}+p.OT{iStep}(:,iMode); 
    end  
    end
    
    %aa=p.T{iStep+1}
   
    
    %combine m with speed restriction
    m{iStep}=min(m{iStep},speedRes);    
    %combine with free driving probability
    freeDriving=ones(length(m{iStep}),1)*freeDrivingProb;
    m{iStep}=combine(m{iStep},freeDriving);
    
    %compute input transitions
    %get nonzero elements of p.T
    ind=find(pTotal.T{iStep+1});
    %obtain input transitions for each state space cell
    for iCell=1:length(ind)
        %compute Gamma matrix
        G=diag(m{iStep}(ind(iCell),:))*Gamma;
        G=normalizeMatrix(G);

        %execute Gamma matrix
        p.T{iStep+1}(ind(iCell),:)=G*p.T{iStep+1}(ind(iCell),:)';
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