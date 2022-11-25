function [mTotal]=combineOptimized(resProb,freeProb,nrOfInputs)
% combines - combines two cell arrays of selection vectors m1 and m2 to a
% total selection vector
%
% Syntax:  
%    [mTotal]=combine(m1,m2)
%
% Inputs:
%    pTotal - total probabilities of different time steps
%    field - partition object
%
% Outputs:
%    m1 - first array of selection vectors
%    m2 - second array of selection vectors
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
% Written:      01-July-2008 
% Last update:  04-November-2008
%               01-July-2009
% Last revision: ---

%------------- BEGIN CODE --------------


%obtain number of states
nrOfStates=length(resProb)/nrOfInputs;

%mode loop
for i=1:nrOfInputs
    %set iMode
    iMode=nrOfInputs-i+1;
    
    %set new state indices
    stateIndices=nrOfInputs*((1:nrOfStates)-1)+iMode;
    
    %find values where the restriction is lower than the free driving
    %probability
    ind=find(resProb(stateIndices)<freeProb(stateIndices));
    ind2=find(resProb(stateIndices)>=freeProb(stateIndices));
    
    %set values for those indices to the restriction probabilities
    selStates=nrOfInputs*(ind-1)+iMode;
    selStates2=nrOfInputs*(ind2-1)+iMode;
    mTotal(selStates)=resProb(selStates);
    mTotal(selStates2)=freeProb(selStates2);
    
    %add difference to free probability
    diff(selStates,1)=freeProb(selStates)-resProb(selStates);
    selStatesNew=nrOfInputs*(ind-1)+iMode-1;
    if iMode>1
        freeProb(selStatesNew)=freeProb(selStatesNew)+diff(selStates);
    end
end

%------------- END OF CODE --------------