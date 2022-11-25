function [posProb,velProb]=project(pTotal,field)
% project - Projects probabilities onto position and velocity probabilities
%
% Syntax:  
%    [posProb,velProb]=project(pTotal,field)
%
% Inputs:
%    pTotal - total probabilities of different time steps
%    field - partition object
%
% Outputs:
%    posProb - position probability
%    velProb - velocity probability
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
% Written:      17-July-2008 
% Last update:  13-November-2017
% Last revision:---

%------------- BEGIN CODE --------------

%obtain number of segments
nrOfSeg=field.nrOfSegments;

%projection matrix for position
P2position=sparse(matrixbuilder(nrOfSeg(1),nrOfSeg(2),0)); 
%projection matrix for velocity
P2velocity=sparse(matrixbuilder(nrOfSeg(1),nrOfSeg(2),1)); 

for iStep=1:length(pTotal)
    %project
    posProb{iStep}=P2position*pTotal{iStep};
    velProb{iStep}=P2velocity*pTotal{iStep};
end

%------------- END OF CODE --------------