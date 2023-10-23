function tp = transitionProbability_reach(R,tranFrac,field)
% transitionProbability_reach - Calculate the transition probability from 
% the actual cell to the reachable cells using reachability analysis.
%
% Syntax:  
%    tp = transitionProbability_reach(niP,tranFrac,field)
%
% Inputs:
%    niP - niP (non intersecting polytopes)
%    tranFrac - transitionFraction (contains probabilities of transitioning 
%       to other mode of the original hybrid automaton)
%    field - partition of the state space 
%
% Outputs:
%    tp - transition probability struct
%
% Example: 
%    -
%
% Other m-files required: tbd
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author:       Matthias Althoff
% Written:      15-September-2006
% Last update:  09-October-2006
%               26-March-2008
%               29-September-2009
%               31-July-2017
%               24-July-2020
% Last revision:---

%------------- BEGIN CODE --------------

%initialize--------------------------------------------------------
nrOfStates = nrOfCells(field);
tp.T(1:(nrOfStates+1),1)=0; %tp: transition probability
tp.OT(1:(nrOfStates+1),1)=0; %tp: transition probability
%------------------------------------------------------------------


%get cells that might intersect with the reachable set-------------
for k=1:length(R)
    % solution for time point
    %polytope conversion if niP{k} is a zonotope
    niP_tp = R(k).timePoint.set{end};
    if isa(niP_tp,'zonotope')
        niP_tp = polytope(niP_tp);
    end
    % intersection probabilities
    [~, iP_tp] = exactIntersectingCells(field, niP_tp);
    % add partial transition probbailities from the considered time
    % interval
    tp.T = tp.T + tranFrac{k}*iP_tp;
    
    % solution for time interval
    % number of sets in one location
    nrOfSets = length(R(k).timeInterval.set);
    for i=1:nrOfSets
        %polytope conversion if niP{k} is a zonotope
        niP_ti = R(k).timeInterval.set{i};
        if isa(niP_ti,'zonotope')
            niP_ti = polytope(niP_ti);
        end
        % intersection probabilities
        [~, iP_ti] = exactIntersectingCells(field, niP_ti);
        % add partial transition probbailities from the considered time
        % interval
        tp.OT = tp.OT + tranFrac{k}/nrOfSets*iP_ti;
    end
end

%------------- END OF CODE --------------