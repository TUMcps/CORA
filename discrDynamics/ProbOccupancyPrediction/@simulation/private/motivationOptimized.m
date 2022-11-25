function [mot]=motivationOptimized(w,sel,p,simOptions)
% motivation - compute the motivation value for driving on a certain lane
%
% Syntax:  
%    [mot]=motivation(p,w)
%
% Inputs:
%    p - probability matrix for the joint probability of states and inputs
%    w - weighting vector
%
% Outputs:
%    mot - motivation value
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
% Written:      15-October-2009
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%get nr of inputs
nrOfInputs=length(simOptions.initialInputProb);

%reshape selection vector
selTmp=reshape(sel,nrOfInputs,prod(size(sel))/nrOfInputs);
sel=selTmp';

%init pNew
pNew=0*sel;

%compute total probability
pTotal=simOptions.projMat*p;
%get nonzero elements of total probability
ind=find(pTotal);
%obtain input transitions for each state space cell
for iCell=1:length(ind)
    %compute new probabilities
    pNew(ind(iCell),:)=sel(ind(iCell),:)*pTotal(ind(iCell));
end

%sum up the probabilities for the different inputs
s=sum(pNew);

%multiply result with the weighting vector
mot=w*s';

%------------- END OF CODE --------------