function [transSet,mergedLabels] = mergeTransitionSets(pHA,locID,varargin)
% mergeTransitionSets - computes the transition set of a combined location 
%    from the transition sets of the subcomponents
%
% Syntax:
%    transSet = mergeTransitionSets(pHA,locID)
%    [transSet,mergedLabels] = mergeTransitionSets(pHA,locID,allLabels)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    locID - IDs of the currently active locations
%    allLabels - information about all synchronization labels
%
% Outputs:
%    transSet - cell array containing the resulting transition sets
%    mergedLabels - list of synchronization labels which were used to
%                   construct a merged transition
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: locationProduct, mergeInvariants

% Authors:       Johann Schoepfer, Niklas Kochdumper, Maximilian Perschl
% Written:       14-June-2018
% Last update:   16-March-2022 (MP, implement synchronization labels)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments... re-compute list of label occurrences if
% necessary (note: we do not check whether the struct is correct...)
narginchk(2,3);
allLabels = setDefaultValues({[]},varargin);
if isnumeric(allLabels) && isempty(allLabels)
    allLabels = labelOccurrences(pHA);
end

% number of components
numComps = length(pHA.components);

% pre-compute list of flows
flowList = arrayfun(@(comp) [comp.location.contDynamics]',...
    pHA.components,'UniformOutput',false);

% compute total input dimension
M = size(vertcat(pHA.bindsInputs{:}),1);
% rewrite input binds analogously to state binds (global counter)
idxCompInput = arrayfun(@(i) size(vertcat(pHA.bindsInputs{1:i-1}),1)+1:...
    size(vertcat(pHA.bindsInputs{1:i}),1),1:numComps,'UniformOutput',false);

% initialize resulting set of transitions
transSet = transition();
% counter for number of transitions (cannot be known in advance due to
% synchronization labels)
cnt = 1;

% list of checked and merged synchronization labels
checkedLabels = {};
mergedLabels = {};

% loop over all subcomponents
for i = 1:numComps

    % read out list of transitions of for currently active location of i-th
    % subcomponent
    transList = pHA.components(i).location(locID(i)).transition;

    % loop over all transitions for the current subcomponent
    for t = 1:length(transList)
        
        % read out transition and synchronization label
        trans = transList(t);
        syncLabel = trans.syncLabel;

        if ismember(syncLabel,checkedLabels)
            % skip labels that have already been checked
            continue;

        elseif isempty(syncLabel)
            % a transition without a synchronization label can be directly
            % converted to a transition of the composed automaton

            % lift to higher dimension, synchronize to resolve input
            % dependencies, and write to list of transitions
            trans_lift = lift(trans,pHA.nrOfStates,M,pHA.bindsStates{i},...
                idxCompInput{i},false);
            % list all states with identity mapping (all but the states
            % of the ith component since it's the only participating
            % component in the resulting transition)
            idStates = setdiff(1:pHA.nrOfStates,pHA.bindsStates{i});
            % read out active flow equation for each component given the
            % currently active location ID vector... required for resolving
            % input binds via the output equations
            flowList_active = arrayfun(@(i) flowList{i}(locID(i)),...
                1:numComps,'UniformOutput',false)';
            transSet(cnt) = synchronize(trans_lift,idStates,locID,i,...
                pHA.bindsStates,pHA.bindsInputs,flowList_active);
        
        else
            % the goal of this branch is to gather all transitions with the
            % same synchronization label, then we call lift on all
            % individual transitions and synchronize them afterward in
            % the same way as above

            % append new synchronization label to list of checked labels
            % to that we can skip further instances later
            checkedLabels = [checkedLabels; syncLabel];
            
            % gather all transitions with the given synchronization
            % label; if not all components which contain the 
            % synchronization label are currently in a location where 
            % such a transition could be taken, the transitions cannot
            % be synchronized (returns empty, goes to 'continue')
            [labelTransSet,labelCompIdx] = ...
                priv_gatherTransitions(pHA,locID,syncLabel,allLabels);
            if isempty(labelTransSet)
                continue
            end

            % save label for product of invariants
            mergedLabels = [mergedLabels; syncLabel];

            % for all non-participating components, we keep their states
            % as they are, and therefore set the respective entry in the
            % vector of all states to 1 (via the state binds)
            idStates = setdiff(1:pHA.nrOfStates,cell2mat(pHA.bindsStates(labelCompIdx)));
            
            % create a new transition from the transitions that are to
            % be synchronized, which models the relevant transitions being
            % executed simultaneously

            % lift all transitions to higher dimension, synchronize, and
            % resolve input dependencies
            for j=1:length(labelTransSet)
                labelTransSet(j) = lift(labelTransSet(j),pHA.nrOfStates,M,...
                    pHA.bindsStates{j},idxCompInput{labelCompIdx(j)},false);
            end
            transSet(cnt) = synchronize(labelTransSet,idStates,locID,...
                labelCompIdx,pHA.bindsStates,pHA.bindsInputs,flowList);
        end

        % increment counter of resulting transitions
        cnt = cnt + 1;
    end
end

% compute derivatives of transition for further use
fpath = [CORAROOT filesep 'models' filesep 'auxiliary' filesep ...
    'parallelHybridAutomaton' filesep pHA.name filesep ...
    'location_' replace(num2str(reshape([1;1;2],1,[])),' ','')];
fname = 'merged';
transSet = derivatives(transSet,fpath,fname);

% ------------------------------ END OF CODE ------------------------------
