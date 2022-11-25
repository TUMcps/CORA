function [pNew,pTotal]=singleStepDriving(simOptions,markovChainSpec,p)
% singleStepDriving - like driving.m, but for a single time step
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
% Written: 30-October-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%load data
m=simOptions.selectionVector;
T=simOptions.transitionMatrix;
Gamma=simOptions.Gamma;
speedRes=simOptions.speedRes;
reachIndices=simOptions.reachIndices;
freeDrivingProb=simOptions.freeDrivingProb;

%obtain the time step
iStep=length(p.T);

%define zero vector z 
z=0*p.T{1}(:,1);

    
%initialize pTotal
pTotal.T=z;
pTotal.OT=z;

%delete probabilities outside the reachable set
ind=reachIndices{iStep}+1; %0 is the outside cell
pBefore=sum(sum(p.T{iStep}));
pTemp=0*p.T{iStep};
pTemp(ind,:)=p.T{iStep}(ind,:);
pTemp(1,:)=p.T{iStep}(1,:); %do not cancel outside probabilities
s=sum(sum(pTemp));
if s~=0
    p.T{iStep}=pTemp/s*pBefore;
end

%input mode loop
for iMode=1:markovChainSpec.nrOfInputs        
    %time step solution
    pNew.T(:,iMode)=T.T{iMode}*p.T{iStep}(:,iMode);
    %time interval solution
    pNew.OT(:,iMode)=T.OT{iMode}*p.T{iStep}(:,iMode);

    %compute total probabilities 
    pTotal.T=pTotal.T+pNew.T(:,iMode);
    pTotal.OT=pTotal.OT+pNew.OT(:,iMode); 
end  

%combine m with speed restriction
m{iStep}=min(m{iStep},speedRes);    
%combine with free driving probability
freeDriving=ones(length(m{iStep}),1)*freeDrivingProb;
m{iStep}=combine(m{iStep},freeDriving);   

%compute input transitions
%get nonzero elements of p.T
ind=find(pTotal.T);
%obtain input transitions for each state space cell
for iCell=1:length(ind)
    %compute Gamma matrix
    G=diag(m{iStep}(ind(iCell),:))*Gamma;
    G=normalizeMatrix(G);

    %execute Gamma matrix
    pNew.T(ind(iCell),:)=G*pNew.T(ind(iCell),:)';
end

    

%------------- END OF CODE --------------