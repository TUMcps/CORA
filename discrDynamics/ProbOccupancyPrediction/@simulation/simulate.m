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

% Authors:       Matthias Althoff
% Written:       17-July-2008 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%compute initial probability
stateField = obj.simOptions.stateField;
initialStateSet = obj.simOptions.initialStateSet;
[~, obj.simOptions.initialProbability] = exactIntersectingCells(stateField, initialStateSet);
%compute Gamma matrices
obj.simOptions.Gamma=priv_gammaMatrix(obj.markovChainSpec,obj.simOptions);
%compute speed restriction
obj.simOptions.speedRes=priv_speedRestriction(obj.simOptions, obj.markovChainSpec);
%compute deterministic reachable cell indices
[obj.simOptions.reachIndices]=priv_reachableSet(obj.simOptions,obj.markovChainSpec);


%pick specialized algorithm for different types of traffic situations
switch obj.simOptions.mode
    case 'autonomousDriving'
        [p,pTotal]=priv_autonomousDriving(obj.simOptions,obj.markovChainSpec);
    case 'freeDriving'
        [p,pTotal]=priv_freeDriving(obj.simOptions,obj.markovChainSpec);
    case 'vehicleFollowing'
        [p,pTotal]=priv_vehicleFollowing(obj.simOptions,obj.markovChainSpec);
    case 'roadCrossing'
        [p,pTotal]=priv_roadCrossing(obj.simOptions,obj.markovChainSpec);
    case 'laneChanging'
        [p,pTotal,lcEvolProb]=priv_laneChanging(obj.simOptions,obj.markovChainSpec);
end
		
if strcmp(obj.simOptions.mode,'laneChanging')
    %project onto position and velocity probabilities
    [posProb.left,velProb.left]=priv_project(pTotal.left.OT,obj.simOptions.stateField); 
    [posProb.right,velProb.right]=priv_project(pTotal.right.OT,obj.simOptions.stateField); 
    %project onto input probabilities
    inputProb.left=priv_inputDist(p.left.OT);
    inputProb.right=priv_inputDist(p.right.OT);
    %store ratio of lane change/lane keeping
    obj.result.lcEvolProb=lcEvolProb;
else
    %project onto position and velocity probabilities
    [posProb,velProb]=priv_project(pTotal.OT,obj.simOptions.stateField);
    [posProb_T,velProb_T]=priv_project(pTotal.T,obj.simOptions.stateField);
    %project onto input probabilities
    inputProb=priv_inputDist(p.OT);    
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

% ------------------------------ END OF CODE ------------------------------
