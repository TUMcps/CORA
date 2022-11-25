function [obj]=simulate(obj)
% simulate - Simulates a traffic participant
%
% Syntax:  
%    [obj]=simulate(obj)
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

% Author: Matthias Althoff
% Written: 17-July-2008 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%compute initial probability
stateField = obj.simOptions.stateField;
initialStateSet = obj.simOptions.initialStateSet;
[~, obj.simOptions.initialProbability] = exactIntersectingCells(stateField, initialStateSet);
%compute Gamma matrices
obj.simOptions.Gamma=gammaMatrix(obj.markovChainSpec,obj.simOptions);
%compute speed restriction
obj.simOptions.speedRes=speedRestriction(obj.simOptions, obj.markovChainSpec);
%compute deterministic reachable cell indices
[obj.simOptions.reachIndices]=reachableSet(obj.simOptions,obj.markovChainSpec);


%pick specialized algorithm for different types of traffic situations
switch obj.simOptions.mode
    case 'autonomousDriving'
        [p,pTotal]=autonomousDriving(obj.simOptions,obj.markovChainSpec);
    case 'freeDriving'
        [p,pTotal]=freeDriving(obj.simOptions,obj.markovChainSpec);
    case 'vehicleFollowing'
        [p,pTotal]=vehicleFollowing(obj.simOptions,obj.markovChainSpec);
    case 'roadCrossing'
        [p,pTotal]=roadCrossing(obj.simOptions,obj.markovChainSpec);
    case 'laneChanging'
        [p,pTotal,lcEvolProb]=laneChanging(obj.simOptions,obj.markovChainSpec);
end
		
if strcmp(obj.simOptions.mode,'laneChanging')
    %project onto position and velocity probabilities
    [posProb.left,velProb.left]=project(pTotal.left.OT,obj.simOptions.stateField); 
    [posProb.right,velProb.right]=project(pTotal.right.OT,obj.simOptions.stateField); 
    %project onto input probabilities
    inputProb.left=inputDist(p.left.OT);
    inputProb.right=inputDist(p.right.OT);
    %store ratio of lane change/lane keeping
    obj.result.lcEvolProb=lcEvolProb;
else
    %project onto position and velocity probabilities
    [posProb,velProb]=project(pTotal.OT,obj.simOptions.stateField);
    [posProb_T,velProb_T]=project(pTotal.T,obj.simOptions.stateField);
    %project onto input probabilities
    inputProb=inputDist(p.OT);    
end

%write results to object
obj.result.p=p;
obj.result.pTotal=pTotal;
obj.result.positionProbability=posProb;
obj.result.velocityProbability=velProb;
obj.result.positionProbability_T=posProb_T;
obj.result.velocityProbability_T=velProb_T;
obj.result.inputProbability=inputProb;

% %plot exemplary result
% field=obj.simOptions.stateField;
% pAll=0*pTotal.T{1};
% for i=1:10
%     pAll=pAll+pTotal.OT{i};
% end
% pAll=pAll/i;
% plotP(field,pAll,'k');

%------------- END OF CODE --------------