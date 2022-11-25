function [p,pTotal]=autonomousDrivingOptimized(simOptions,markovChainSpec)
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
% Last update:  13-November-2017
% Last revision:---

%------------- BEGIN CODE --------------

%get field
field=simOptions.stateField;

%position, velocity and mode vector
posVec=simOptions.autonomousCarTrajectory.position;
velVec=simOptions.autonomousCarTrajectory.velocity;
modeVec=simOptions.autonomousCarTrajectory.mode;

% map results to uniform probabilities
% assume that car is x meters long
% assume that speed is uncertain within y m/s
if isfield(simOptions,'uncertainEgoPos')
    uncertainPos=simOptions.uncertainEgoPos;
else
    uncertainPos=3;
end
if isfield(simOptions,'uncertainEgoVel')
    uncertainVel=simOptions.uncertainEgoVel;
else
    uncertainVel=1;
end

for iStep=1:simOptions.runs
    
    %generate interval hull for time point solution
    IHtp=interval([posVec(iStep)-uncertainPos; velVec(iStep)-uncertainVel],...
        [posVec(iStep)+uncertainPos; velVec(iStep)+uncertainVel]);
    
    %get interval hull for time interval
    minPos=min(posVec(iStep),posVec(iStep+1));
    maxPos=max(posVec(iStep),posVec(iStep+1));
    
    minVel=min(velVec(iStep),velVec(iStep+1));
    maxVel=max(velVec(iStep),velVec(iStep+1));    
    
    IHti=interval([minPos-uncertainPos; minVel-uncertainVel], ...
        [maxPos+uncertainPos; maxVel+uncertainVel]);
    
    %convert intervals to probabilities
    %time point
    [~, pTP] = exactIntersectingCells(field,IHtp);
    %time interval
    [~, pTI] = exactIntersectingCells(field,IHti);
    
    %total probability = partial probability
    pTotal.T{iStep}=pTP;
    pTotal.OT{iStep}=pTI;
    
    % mode dependent probabilities
    for iMode=1:markovChainSpec.nrOfInputs
        if iMode==modeVec(iStep)
            p.T{iStep}(:,iMode)=pTP;
            p.OT{iStep}(:,iMode)=pTI;             
        else
            p.T{iStep}(:,iMode)=0*pTP;
            p.OT{iStep}(:,iMode)=0*pTI;  
        end
    end
    
    %reshape probability distribution
    pTmp.T=reshape(p.T{iStep}',prod(size(p.T{iStep})),1);
    pTmp.OT=reshape(p.OT{iStep}',prod(size(p.OT{iStep})),1);
    p.T{iStep}=[];
    p.OT{iStep}=[];
    p.T{iStep}=pTmp.T; 
    p.OT{iStep}=pTmp.OT; 
end


%------------- END OF CODE --------------