function [mTotal]=combine(resProb,freeProb)
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

% Author: Matthias Althoff
% Written: 01-July-2008 
% Last update: 04-November-2008
% Last revision: ---

%------------- BEGIN CODE --------------

%continue here!!

%obtain number of inputs
nrOfInputs=length(resProb(1,:));

%mode loop
for i=1:nrOfInputs
    %set iMode
    iMode=nrOfInputs-i+1;
    
    %find values where the restriction is lower than the free driving
    %probability
    ind=find(resProb(:,iMode)<freeProb(:,iMode));
    ind2=find(resProb(:,iMode)>=freeProb(:,iMode));
    
    %set values for those indices to the restriction probabilities
    mTotal(ind,iMode)=resProb(ind,iMode);
    mTotal(ind2,iMode)=freeProb(ind2,iMode);
    
    %add difference to free probability
    diff(ind,1)=freeProb(ind,iMode)-resProb(ind,iMode);
    if iMode>1
        freeProb(ind,iMode-1)=freeProb(ind,iMode-1)+diff(ind,1);
    end
end

%------------- END OF CODE --------------