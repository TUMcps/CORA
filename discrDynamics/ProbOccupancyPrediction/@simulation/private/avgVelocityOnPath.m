function [avgVelocity]=avgVelocityOnPath(pTotal,field)
% avgVelocityOnPath - computes the average velocity on path segments
%
% Syntax:  
%    [avgVelocity]=avgVelocityOnPath(pTotal,field)
%
% Inputs:
%    pTotal - total probabilities of different time steps
%    field - partition object
%
% Outputs:
%    avgVelocity - average velocity on path
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
% Last update:  13-November-2017
% Last revision:---

%------------- BEGIN CODE --------------

%obtain number of segments
nrOfSeg=field.nrOfSegments;
segLength=(field.intervals(:,2)-field.intervals(:,1))./nrOfSeg;

%projection matrix for position
P2position=sparse(matrixbuilder(nrOfSeg(1),nrOfSeg(2),0)); 

%compute center velocities for each velocity interval
centerVel=0.5*segLength(2):segLength(2):nrOfSeg(2)*segLength(2)-0.5*segLength(2);

%weight projection matrix with center velocities
P2avgVel(:,1)=P2position(:,1);
for iPos=1:nrOfSeg(1)
    for iVel=1:nrOfSeg(2)
        col=1+(iVel-1)*nrOfSeg(1)+iPos;
        P2avgVel(:,col)=centerVel(iVel)*P2position(:,col);
    end
end

%compute average total probability
pTotalSum=0*pTotal{1};
for iStep=1:length(pTotal)
    %sum up
    pTotalSum=pTotalSum+pTotal{iStep};
end

%compute average probability
pTotalAvg=pTotalSum/sum(pTotalSum);

%obtain average velocity
avgVelocity=P2avgVel*pTotalAvg;
pathProb=P2position*pTotalAvg;
nonZeroInd=find(pathProb);
avgVelocity(nonZeroInd)=avgVelocity(nonZeroInd)./pathProb(nonZeroInd);

%------------- END OF CODE --------------