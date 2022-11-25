function [inputProb]=inputDistOptimized(p,projMat)
% inputDistOptimized - Returns the input probability distribution of the simulated
% traffic participant
%
% Syntax:  
%    [inputProb]=inputDist(p)
%
% Inputs:
%    p - joint probabilities (state and input) of different time steps
%
% Outputs:
%    inputProb - probability distribution of the input
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      09-October-2009
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%obtain number of inputs and states
nrOfStates=length(projMat(:,1));
nrOfInputs=length(projMat(1,:))/nrOfStates;

%obtain projection matrix
projMat=matrixbuilder(nrOfInputs,nrOfStates,0);
projMat(:,1)=[];

for iStep=1:length(p)
    %project
    inputProb{iStep}=projMat*p{iStep};
end

%------------- END OF CODE --------------