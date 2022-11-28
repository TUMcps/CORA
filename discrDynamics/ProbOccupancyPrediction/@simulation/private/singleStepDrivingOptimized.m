function [pNew,pNewTotal]=singleStepDrivingOptimized(simOptions,p)
% singleStepDrivingOptimized - like drivingOptimized.m, but for a single 
% time step
%
% Syntax:  
%    [pNew,pNewTotal]=singleStepDrivingOptimized(simOptions,p)
%
% Inputs:
%    simOptions - simulation options
%    p - cell array of probability vectors
%
% Outputs:
%    pNew - new probability distributions for next time point
%    pNewTotal - new total probability distributions for next time point
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
mFull=simOptions.selectionVector;
T=simOptions.transitionMatrix;
projMat=simOptions.projMat;
GammaFull=simOptions.GammaFull;
speedRes=simOptions.speedRes;
freeDrivingProb=simOptions.freeDrivingProb;

%nr of inputs and states
nrOfInputs=length(freeDrivingProb);
nrOfStates=length(simOptions.initialProbability);
nrOfCombinedStates=length(p.T{1});

%reshape free driving probability distribution
freeDriving=ones(nrOfStates,1)*freeDrivingProb;
freeDrivingTmp=reshape(freeDriving',prod(size(freeDriving)),1);
freeDriving=[];
freeDriving=freeDrivingTmp;

%obtain the time step
iStep=length(p.T);

%time step solution
p.T{iStep+1}=T.T*p.T{iStep};
%time interval solution
p.OT{iStep}=T.OT*p.T{iStep};

%compute total probabilities 
pTotal.T{iStep+1}=projMat*p.T{iStep+1};
pTotal.OT{iStep}=projMat*p.OT{iStep};  


%get nonzero indices of the reachable set
nonZeroIndTmp=find(pTotal.T{iStep+1});
if length(nonZeroIndTmp)>0
    for i=1:length(nonZeroIndTmp)
        nonZeroIndTotal(((i-1)*nrOfInputs+1):(i*nrOfInputs))=(nonZeroIndTmp(i)-1)*nrOfInputs+(1:nrOfInputs);
    end

    %initialize mAct to 0
    mPartial=sparse(nrOfCombinedStates,1);    
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

%write results to pNew and pNewTotal
pNew.T=p.T{iStep+1};
pNew.OT=p.OT{iStep};
pNewTotal.T=pTotal.T{iStep+1};
pNewTotal.OT=pTotal.OT{iStep};


%------------- END OF CODE --------------