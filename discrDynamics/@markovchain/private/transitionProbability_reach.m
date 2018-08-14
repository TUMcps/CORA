function tp = transitionProbability_reach(niP,tranFrac,field)
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
% Last revision:---

%------------- BEGIN CODE --------------

%initialize--------------------------------------------------------
nrOfStates = nrOfCells(field);
tp(1:(nrOfStates+1),1)=0; %tp: transition probability
%------------------------------------------------------------------


%get cells that might intersect with the reachable set-------------
for k=1:length(niP)
    for i=1:length(niP{k})
        if ~iscell(niP{k}{i})
            %polytope conversion if niP{k}{i} is a zonotope
            if isa(niP{k}{i},'zonotope')
                niP{k}{i}=polytope(niP{k}{i});
            end
            % intersection probabilities
            [~, iP] = exactIntersectingCells(field,niP{k}{i});
            % add partial transition probbailities from the considered time
            % interval
            tp=tp+tranFrac{k}/length(niP{k})*iP;
        else
            for j=1:length(niP{k}{i})
                %polytope conversion if niP{k}{i} is a zonotope
                if isa(niP{k}{i}{j},'zonotope') || isa(niP{k}{i}{j},'zonotopeBundle')
                    niP{k}{i}{j}=polytope(niP{k}{i}{j});
                    %intersection
                    [~, iP] = exactIntersectingCells(field,niP{k}{i}{j});
                elseif isa(niP{k}{i}{j}.set,'zonotope') || isa(niP{k}{i}{j}.set,'zonotopeBundle')
                    niP{k}{i}{j}.set=polytope(niP{k}{i}{j}.set);
                    %intersection
                    [~, iP] = exactIntersectingCells(field,niP{k}{i}{j}.set);
                end

                % add transition probabilities
                tp=tp+tranFrac{k}/length(niP{k})*iP;
            end
        end
    end
end

%------------- END OF CODE --------------