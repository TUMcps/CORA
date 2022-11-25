function [T, projMat, GammaFull] = convertTransitionMatrix(obj, gamma)
% convertTransitionMatrix - converts the transition matrix of a
% Markov-chain such that it can be used for an optimized update as
% presented in M. Althoff, O. Stursberg, and M. Buss. Model-based 
% probabilistic collision detection in autonomous driving. 
% IEEE Transactions on Intelligent Transportation Systems, 10:299 - 310, 2009.
%
% Syntax:  
%    [T, projMat, GammaFull] = convertTransitionMatrix(obj)
%
% Inputs:
%    obj - markovchain object (unoptimized structure)
%
% Outputs:
%    T - transition matrix (optimized structure)
%    projMat - projection matrix
%    GammaFull - full Gamma matrix (see 2009 Tansactions on Intelligent Transportation Systems paper)
%    
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
%               01-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%retrieve transition matrix from Markov Chain
T=get(obj,'T');

%compute projection matrix
disp('compute projMaT');
projMat=projectionMatrix(T.T);

%compute full Gamma-matrix
disp('compute Gamma');
GammaFull=computeGammaFull(T.T, gamma);

%convert
disp('convert T');
T.T=convert(T.T);
disp('convert OT');
T.OT=convert(T.OT);


function T=convert(Phi)

%obtain number of inputs and states
nrOfStates=length(Phi{1});
nrOfInputs=length(Phi);

%allocate sparse Matrix and reserve space for nnz([Phi{:}]) entries
T=spalloc(nrOfInputs*nrOfStates,nrOfInputs*nrOfStates,nnz([Phi{:}]));

%reorganize structure
 for iInput=1:nrOfInputs
     disp(['Input ',num2str(iInput)]);
      
     %get nonzero entries
    [goalInd,initialInd]=find(Phi{iInput});
 
    %New row in the T matrix
    i=iInput+nrOfInputs*(goalInd-1);
    
    %New collum in the T matrix
    j=iInput+nrOfInputs*(initialInd-1);
    
    %The entry that shall be copied to T(i,j) (as vector)
    s=Phi{iInput}(sub2ind(size(Phi{iInput}),goalInd,initialInd));
    
    T=T+sparse(i,j,s,nrOfInputs*nrOfStates,nrOfInputs*nrOfStates);
 end


function GammaFull=computeGammaFull(Phi, gamma)

%obtain number of inputs and states
nrOfStates=length(Phi{1});
nrOfInputs=length(Phi);
nrOfCombinedStates=nrOfInputs*nrOfStates;

%build Gamma matrix
for i=1:nrOfInputs
    for j=1:nrOfInputs
        Gamma(i,j)=1/((i-j)^2+gamma);
    end
end
Gamma=normalizeMatrix(Gamma);

%obtain GammaFull
GammaFull=sparse(nrOfCombinedStates,nrOfCombinedStates);
for index=1:nrOfStates
    fullIndices=(1:nrOfInputs)+nrOfInputs*(index-1);
    GammaFull(fullIndices,fullIndices)=Gamma;
end


function projMat=projectionMatrix(Phi)

%obtain number of inputs and states
nrOfStates=length(Phi{1});
nrOfInputs=length(Phi);

%obtain projection matrix
projMat=matrixbuilder(nrOfInputs,nrOfStates,1);
projMat(:,1)=[];



%------------- END OF CODE --------------