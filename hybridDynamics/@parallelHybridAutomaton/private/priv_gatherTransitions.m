function [transSet,compIdx,targets] = priv_gatherTransitions(pHA,locID,syncLabel,allLabels)
% priv_gatherTransitions - gather all transitions of a parallel hybrid
%    automaton with the same synchronization label
%
% Syntax:
%    [transSet,compIdx,targets] = priv_gatherTransitions(pHA,locID,syncLabel,allLabels)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    locID - IDs of the currently active locations
%    syncLabel - synchronization label
%    allLabels - information about all synchronization labels
%
% Outputs:
%    transSet - class array of transitions with same synchronization label
%    compIdx - index of components in parallel hybrid automaton where these
%              transition come from
%    targets - target location after transition in each subcomponent
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mergeTransitionSets

% Authors:       Maximilian Perschl, Mark Wetzlinger
% Written:       17-February-2022
% Last update:   ---
% Last revision: 10-October-2024 (MW, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% initialize list of transitions with same label and the index of the
% component they come from, as well as the target location within the
% respective component
transSet = []; compIdx = []; targets = [];

% loop over all components
for i=1:length(pHA.components)

    % read out transitions for currently active location
    transList = pHA.components(i).location(locID(i)).transition;

    % find given synchronization label in the set of transitions
    idx = arrayfun(@(trans) strcmp(trans.syncLabel,syncLabel),...
        transList,'UniformOutput',true);

    % note that, by construction -- checked in the constructor of the
    % parallelHybridAutomaton object -- only one transition per component
    % may have a given synchronization label, so idx is either empty or of
    % length 1 -> we can use this index to read out the relevant transition
    if ~isempty(idx) && any(idx)
        transSet = [transSet; transList(idx)];
        compIdx = [compIdx; i];
        targets = [targets; transList(idx).target];
    end
end

% the transition is only active if in all components where the label
% occurs, one of these transitions is active
% (the lengthy variable on the right-hand side computes the number of
% components in which <label> occurs)
if length(transSet) ~= ...
        length(unique(allLabels(strcmp({allLabels.name},syncLabel)).component))
    transSet = []; compIdx = []; targets = [];
end

% ------------------------------ END OF CODE ------------------------------
