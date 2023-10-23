function [mergedTrans,tracker] = instantTransition(pHA,locID,allLabels,tracker)
% instantTransition - if the currently active locations allow for an
%    instant transition, i.e.
%       - either there is a transition without a synchronization label, or
%       - all transitions with the same synchronization label can be taken
%         instantly,
%    this transition is constructed here; any subsequent instant
%    transitions are executed in the next iteration of the main loop
%
% Syntax:
%    [mergedTrans,tracker] = instantTransition(pHA,list,locID,allLabels,tracker)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    locID - vector of location IDs for all components
%    allLabels - struct describing synchronization labels
%    tracker (struct) - meta data about the current run
%
% Outputs:
%    mergedTrans - merged transition
%    tracker (struct) - meta data about the current run
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: parallelHybridAutomaton/reach, labelOccurrences

% Authors:       Mark Wetzlinger
% Written:       20-June-2022
% Last update:   ---
% Last revision: 06-June-2023 (MW, return merged transition)

% ------------------------------ BEGIN CODE -------------------------------

% logical indices for instant transitions in currently active locations
immTransNoLabel = cell(length(pHA.components),1);
immTransLabel = cell(length(pHA.components),1);

% loop over all components, check currently active location
for i=1:length(pHA.components)

    % transition of location
    transList = pHA.components(i).location(locID(i)).transition;
    immTransNoLabel{i} = false(length(transList),1);
    immTransLabel{i} = false(length(transList),1);

    % loop over transitions in current location, check for fullspace guards
    % with or without synchronization labels
    for j=1:length(transList)
        if isa(transList(j).guard,'fullspace')
            if isempty(transList(j).syncLabel)
                immTransNoLabel{i}(j) = true;
            else
                immTransLabel{i}(j) = true;
            end
        end
    end

end

% order for instant transitions:
% 1. if there are instant transitions without synchronization labels,
%    these are transitions are merged and taken first
% 2. if there are no instant transitions without synchronization labels,
%    we have to go over all instant transitions with synchronization
%    labels by checking for each label if there is an active transition in
%    each component where the label occurs

if ~any(cell2mat([immTransNoLabel; immTransLabel]))
    % no instant transitions at all

    % since only instant transitions could cause a livelock and the
    % tracker is currently only used for that, we can skip the entries for
    % the transition and synchronization labels this time (use dummy
    % values which cannot yield a livelock)
    tracker(end).transition = zeros(1,length(locID));
    tracker(end).syncLabel = '';

    % quick exit
    mergedTrans = {}; return

elseif any(cell2mat(immTransNoLabel))
    % there are instant transitions without synchronization labels
    % -> can always be constructed and taken

    % update synchronization label in tracker before init of dummy label
    tracker(end).syncLabel = '';

    % define dummy label so that merging of transitions results in
    % one transition, extend list of all labels (will not affect
    % the general list as this is not a return value)
    syncLabel = 'unifiedTransWithoutSyncLabels';
    allLabels = [allLabels; struct('name',syncLabel,'component',[],...
        'location',[],'transition',[],'instant',[])];

    % check all current locations for instant transitions...
    % note: by definition, there can only be one outgoing instant
    % transition without a synchronization label per location!
    for i=1:length(immTransNoLabel)
        if any(immTransNoLabel{i})
            allLabels(end).component = [allLabels(end).component; i];
            allLabels(end).location = [allLabels(end).location; locID(i)];
            allLabels(end).transition = [allLabels(end).transition; ...
                find(immTransNoLabel{i},1,'first')];
            allLabels(end).instant = [allLabels(end).instant; true];
        end
    end

else % here: any(cell2mat(immTransLabel)) = true

    % assume no synchronization label meets necessary criteria
    syncTransitionFound = false;

    % go over instant transitions with synchronization labels, check if
    % they can be executed -> if yes, find all transition for merge
    for i=1:length(immTransLabel)
        for j=1:length(immTransLabel{i})
    
            % only check instant transitions
            if immTransLabel{i}(j)
    
                % synchronization label of instant transition
                syncLabel = pHA.components(i).location(locID(i)).transition(j).syncLabel;
        
                % skip if synchronization label does not meet criteria
                if ~aux_skipTransition(i,locID,syncLabel,allLabels)
                    syncTransitionFound = true;
                    break
                end

            end
        end

        % take first active synchronization label... if others are active
        % too, they have to wait until the next iteration in the main loop
        if syncTransitionFound
            break
        end
    end

    if ~syncTransitionFound
        % no instant transition with synchronization labels found

        % dummy values for tracker (see if-branch: no instant transitions)
        tracker(end).transition = zeros(1,length(locID));
        tracker(end).syncLabel = '';

        % quick exit
        mergedTrans = {}; return
    end

    % update tracker: used synchronization label
    tracker(end).syncLabel = syncLabel;

end

% from here on, we will construct a transition and that transition
% will then be taken

% logical idx for synchronization label in list of all labels
idxSyncLabel = strcmp({allLabels(:).name},syncLabel);

% initialize components which require identity resets: true (other
% components take transition corresponding to synchronization label)
idxIdReset = true(length(pHA.components),1);

% sychronization label given: already checked that all other
% transitions with same synchronization label are also
% instant transitions; now get these transitions
for k=1:length(pHA.components)
    if any(allLabels(idxSyncLabel).component == k)
        idxIdReset(k) = false;
    end
end

% initialize transition for full state
transList = cell(length(pHA.components),1);

% loop over all components for full transition
for k=1:length(pHA.components)

    % check for identity reset or not
    if idxIdReset(k)
        % append virtual self-transition with fullspace guard set and
        % identity reset function
        stateDim = pHA.components(k).location(locID(k)).contDynamics.dim;
        resetStruct.A = eye(stateDim);
        resetStruct.c = zeros(stateDim,1);
        transList{k,1} = transition(fullspace(stateDim),resetStruct,locID(k),syncLabel);

        % append to list of all labels
        allLabels(idxSyncLabel).component = [allLabels(idxSyncLabel).component; k];
        allLabels(idxSyncLabel).location = [allLabels(idxSyncLabel).location; locID(k)];
        % index this virtual transition by "0" (consistent with livelock check)
        allLabels(idxSyncLabel).transition = [allLabels(idxSyncLabel).transition; 0];
        allLabels(idxSyncLabel).instant = [allLabels(idxSyncLabel).instant; true];

        % update tracker: idx of transition in transition list of location
        tracker(end).transition = [tracker(end).transition; 0];
        
    else
        
        % read out index of transition in list of transition of the
        % k-th component
        idx = allLabels(idxSyncLabel).transition(...
                allLabels(idxSyncLabel).component == k ...
                & allLabels(idxSyncLabel).location == locID(k));

        % take transition with same synchronization label, but
        % remove label for call of mergeTransitionSets below
        transList{k,1} = transition(...
            pHA.components(k).location(locID(k)).transition(idx).guard,...
            pHA.components(k).location(locID(k)).transition(idx).reset,...
            pHA.components(k).location(locID(k)).transition(idx).target,...
            syncLabel);

        % update tracker: idx of transition in transition list of location
        tracker(end).transition = [tracker(end).transition; idx];

    end
    
end

% has the transition-to-be-merged been computed before?
if ~isempty(pHA.mergedTrans)
    % check whether
    % 1) there has been a merged instant transition starting from the same
    %    locID as currently 
    % 2) the same transitions from the transition list have been merged
    sameTrans = all([pHA.mergedTrans.locID] == locID,1) ...
        & all([pHA.mergedTrans.transID] == tracker(end).transition,1);
    if any(sameTrans)
        mergedTrans = pHA.mergedTrans(sameTrans).transition;
        return
    end
end

% construct transition for full state (note that the synchronization label
% is removed during the process)
mergedTrans = mergeTransitionSets(pHA,transList,locID,allLabels);

end


% Auxiliary functions -----------------------------------------------------

function skip = aux_skipTransition(currComp,locID,syncLabel,allLabels)
% check whether current transition is an instant transition: we have to
% go over all other instances of the synchronization label in other
% components to ensure that at least one instance per component is active
% at the moment (i.e., we are in a location with a transition with that 
% synchronization label)

    % logical idx for synchronization label in list of all labels
    idxSyncLabel = strcmp({allLabels(:).name},syncLabel);

    % assume all checks to be ok
    skip = false;

    % components which contain locations with transition with the
    % current synchronization label
    compWithSameLabel = unique(allLabels(idxSyncLabel).component);

    % are all those components in a location with that label and
    % does that transition have an empty guard set?
    for k=1:length(compWithSameLabel)
        % read out component number
        comp = compWithSameLabel(k);

        % skip if number is same as current component
        if comp == currComp
            continue;
        end

        % all entries for given component
        compIdx = allLabels(idxSyncLabel).component == comp;

        % check if any of the corresponding locations is currently active
        if ~any(allLabels(idxSyncLabel).location(compIdx) == locID(comp))
            skip = true; break;
        end

        % check if the corresponding transition is instant
        locIdx = compIdx & ...
            allLabels(idxSyncLabel).location == locID(comp);
        if isempty(locIdx) || ~allLabels(idxSyncLabel).instant(locIdx)
            skip = true; break;
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
