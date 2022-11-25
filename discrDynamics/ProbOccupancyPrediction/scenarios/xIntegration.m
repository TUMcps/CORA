function [pX]=xIntegration(pos,xSegment,intInterval)
% xIntegration - integrates the probability that a car with probability
% distribution p crosses an intersection at the position segment xSegment
% for a certain number of time steps
%
% Syntax:  
%    [pX]=xIntegration(p,xSegment,timeSteps,field)
%
% Inputs:
%    pos - probability distribution of the position of the vehicle
%    xSegment - position segment of the crossing
%    intInterval - integration interval
%
% Outputs:
%    pX - probability that the car has crossed within the specified time
%    Steps 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      19-November-2007 
% Last update:  18-June-2008
% Last revision:---

%------------- BEGIN CODE --------------

for iStep=1:length(pos)
    %generate step interval
    stepInterval(1)=iStep;
    stepInterval(2)=iStep+intInterval;
    %check step interval
    if stepInterval(2)>length(pos)
        stepInterval(2)=length(pos);
    end

    %position probabilities at first and last time step
    pFirst=pos{stepInterval(1)};
    pLast=pos{stepInterval(2)};

    %sum up probabilities in front of the crossing
    pSumFirst=sum(pFirst(1:xSegment));
    pSumLast=sum(pLast(1:xSegment));

    %calculate pX by substracting pSumLast from pSumFirst
    pX(iStep)=pSumFirst-pSumLast;
end

%------------- END OF CODE --------------