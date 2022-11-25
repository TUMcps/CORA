function [pPos]=autonomousDrivingOptimizedPosOnly(simOptions,markovChainSpec)
% autonomousDriving - simulates the autonomous car, based on the 
% autonomous car trajectory
%
% Syntax:  
%    [p,pTotal]=autonomousDriving(simOptions,markovChainSpec)
%
% Inputs:
%    simOptions - simulation options
%    markovChainSpec - Markov-Chain specifications
%
% Outputs:
%    p - probability distributions for different inputs, times
%    pTotal - total probabilities of different times
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      14-October-2009
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%get field
field=simOptions.stateField;

%create new field with pos only
interval=get(field,'intervals');
nrOfSegments=get(field,'nrOfSegments');
field=partition(interval(1,:),nrOfSegments(1));

%position, velocity and mode vector
posVec=simOptions.autonomousCarTrajectory.position;

% map results to uniform probabilities
% assume that car is x meters long
% assume that speed is uncertain within y m/s
if isfield(simOptions,'uncertainEgoPos')
    uncertainPos=simOptions.uncertainEgoPos;
else
    uncertainPos=3;
end


for iStep=1:simOptions.runs
    
    %generate interval hull for time point solution
    IHtp=interval(posVec(iStep)-uncertainPos,posVec(iStep)+uncertainPos);
    
    %get interval hull for time interval
    minPos=min(posVec(iStep),posVec(iStep+1));
    maxPos=max(posVec(iStep),posVec(iStep+1));   
    
    IHti=interval(minPos-uncertainPos,maxPos+uncertainPos);
    
    %convert intervals to probabilities
    %time point
    pPos.T{iStep}=cellIntersection3(field,IHtp);
    %time interval
    pPos.OT{iStep}=cellIntersection3(field,IHti);
end


%------------- END OF CODE --------------