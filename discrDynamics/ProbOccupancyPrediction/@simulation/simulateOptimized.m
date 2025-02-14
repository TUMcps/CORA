function [obj]=simulateOptimized(obj)
% simulateOptimized - Simulates a traffic participant
%
% Syntax:
%    [obj]=simulateOptimized(obj)
%
% Inputs:
%    obj - simulation object
%
% Outputs:
%    obj - simulation object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       17-July-2008 
% Last update:   04-November-2009
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%compute initial probability
stateField = obj.simOptions.stateField;
initialStateSet = obj.simOptions.initialStateSet;
[~, obj.simOptions.initialProbability] = exactIntersectingCells(stateField, initialStateSet);

%compute speed restriction
if isfield(obj.simOptions,'speedRes') && numel(obj.simOptions.speedRes)>0
    disp('Warning: using precalculated speedRes');
else
    speedRes=priv_speedRestrictionOptimized(obj.simOptions, obj.markovChainSpec);
    obj.simOptions.speedRes=speedRes;
end
%compute deterministic reachable cell indices
%[obj.simOptions.reachIndices,minPos,maxPos]=priv_reachableSet(obj.simOptions,obj.markovChainSpec);


%pick specialized algorithm for different types of traffic situations
switch obj.simOptions.mode
    case 'autonomousDriving'
        [p,pTotal]=priv_autonomousDrivingOptimized(obj.simOptions,obj.markovChainSpec);
    case 'autonomousDrivingPosOnly'
        [pProb]=priv_autonomousDrivingOptimizedPosOnly(obj.simOptions,obj.markovChainSpec);
    case 'freeDriving'
        [p,pTotal]=priv_freeDrivingOptimized(obj.simOptions,obj.markovChainSpec);
    case 'vehicleFollowing'
        [p,pTotal]=priv_vehicleFollowingOptimized(obj.simOptions,obj.markovChainSpec);
    case 'roadCrossing'
        [p,pTotal]=priv_roadCrossingOptimized(obj.simOptions,obj.markovChainSpec);
    case 'laneChanging'
        [p,pTotal,lcEvolProb]=priv_laneChangingOptimized(obj.simOptions,obj.markovChainSpec);
end

		
if strcmp(obj.simOptions.mode,'laneChanging')
    %project onto position and velocity probabilities
    [posProb.left,velProb.left]=priv_project(pTotal.left.OT,obj.simOptions.stateField); 
    [posProb.right,velProb.right]=priv_project(pTotal.right.OT,obj.simOptions.stateField); 
    %project onto input probabilities
    inputProb.left=priv_inputDistOptimized(p.left.OT,obj.simOptions.projMat);
    inputProb.right=priv_inputDistOptimized(p.right.OT,obj.simOptions.projMat);
    %obtain average velocity
    avgVel.left=priv_avgVelocityOnPath(pTotal.left.OT,obj.simOptions.stateField);
    avgVel.right=priv_avgVelocityOnPath(pTotal.right.OT,obj.simOptions.stateField);
    %store ratio of lane change/lane keeping
    obj.result.lcEvolProb=lcEvolProb; 
    %dummy for posProb_T, velProb_T
    posProb_T = [];
    velProb_T = [];
elseif strcmp(obj.simOptions.mode,'autonomousDrivingPosOnly')
    p=[]; pTotal=[]; velProb=[]; velProb_T=[]; inputProb=[]; avgVel=[];
    posProb=pProb.OT; posProb_T=pProb.T;
else
    %project onto position and velocity probabilities
    [posProb,velProb]=priv_project(pTotal.OT,obj.simOptions.stateField);
    [posProb_T,velProb_T]=priv_project(pTotal.T,obj.simOptions.stateField);
    %project onto input probabilities
    inputProb=priv_inputDistOptimized(p.OT,obj.simOptions.projMat); 
    %obtain average velocity
    avgVel=priv_avgVelocityOnPath(pTotal.OT,obj.simOptions.stateField); 
end

%write results to object
obj.result.p=p;
obj.result.pTotal=pTotal;

obj.result.positionProbability=posProb;
obj.result.velocityProbability=velProb;
obj.result.positionProbability_T=posProb_T;
obj.result.velocityProbability_T=velProb_T;

obj.result.inputProbability=inputProb;
obj.result.avgVelocity=avgVel;

% %plot exemplary result
% field=obj.simOptions.stateField;
% pAll=0*pTotal.T{1};
% for i=1:10
%     pAll=pAll+pTotal.OT{i};
% end
% pAll=pAll/i;
% plotP(field,pAll,'k');

% ------------------------------ END OF CODE ------------------------------
