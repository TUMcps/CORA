function [list,tracker,restart] = immediateTransition(pHA,list,locID,...
    allLabels,fullR0,fullU,tracker,verbose)
% immediateTransition - if the currently active locations allow for an
%    immediate transition, i.e.
%       - either there is a transition without a synchronization label, or
%       - all transitions with the same synchronization label can be taken
%         immediately,
%    this transition is constructed and executed here; any subsequent
%    immediate transitions are executed in the next iteration of the main
%    loop
%
%    The first entry of 'list' is altered so that the immediate transition
%    does not cause any change in the stored reachable set; the flag
%    'restart' describes whether an immediate transition has been taken.
%
% Syntax:  
%    [list,tracker,restart] = immediateTransition(pHA,list,locID,...
%                               allLabels,fullR0,fullU,tracker,verbose)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    list - queue for reachable set computation
%    locID - vector of location IDs for all components
%    allLabels - struct describing synchronization labels
%    fullR0 - start set of current iteration (full automaton)
%    fullU - input set of current iteration (full automaton)
%    tracker (struct) - meta data about the current run
%    verbose - true/false for logging information on command window
%
% Outputs:
%    list - updated queue for reachable set computation
%    tracker (struct) - meta data about the current run
%    restart - true/false whether immediate transition has been taken
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: parallelHybridAutomaton/reach, labelOccurrences

% Author:       Mark Wetzlinger
% Written:      20-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% logical indices for immediate transitions in currently active locations
immTransNoLabel = cell(length(pHA.components),1);
immTransLabel = cell(length(pHA.components),1);

% loop over all components, check currently active location
for i=1:length(pHA.components)

    % transition of location
    transList = pHA.components(i).location(locID(i)).transition;
    immTransNoLabel{i} = false(length(transList),1);
    immTransLabel{i} = false(length(transList),1);

    % loop over transitions in current location, check for empty guard sets
    % with or without synchronization labels
    for j=1:length(transList)
        if isnumeric(transList(j).guard) && isempty(transList(j).guard)
            if isempty(transList(j).syncLabel)
                immTransNoLabel{i}(j) = true;
            else
                immTransLabel{i}(j) = true;
            end
        end
    end

end

% order for immediate transitions:
% 1. if there are immediate transitions without synchronization labels,
%    these are transitions are merged and taken first
% 2. if there are no immediate transitions without synchronization labels,
%    we have to go over all immediate transitions with synchronization
%    labels by checking for each label if there is an active transition in
%    each component where the label occurs

if ~any(cell2mat([immTransNoLabel; immTransLabel]))
    % no immediate transitions at all

    % since only immediate transitions could cause a livelock and the
    % tracker is currently only used for that, we can skip the entries for
    % the transition and synchronization labels this time (use dummy
    % values which cannot yield a livelock)
    tracker(end).transition = zeros(1,length(locID));
    tracker(end).syncLabel = '';

    % quick exit
    restart = false; return

elseif any(cell2mat(immTransNoLabel))
    % there are immediate transitions without synchronization labels
    % -> can always be constructed and taken

    % update synchronization label in tracker before init of dummy label
    tracker(end).syncLabel = '';

    % define dummy label so that merging of transitions results in
    % one transition, extend list of all labels (will not affect
    % the general list as this is not a return value)
    syncLabel = 'unifiedTransWithoutSyncLabels';
    allLabels = [allLabels; struct('name',syncLabel,'component',[],...
        'location',[],'transition',[],'guardempty',[])];

    % check all current locations for immediate transitions...
    % note: by definition, there can only be one outgoing immediate
    % transition without a synchronization label per location!
    for i=1:length(immTransNoLabel)
        if any(immTransNoLabel{i})
            allLabels(end).component = [allLabels(end).component; i];
            allLabels(end).location = [allLabels(end).location; locID(i)];
            allLabels(end).transition = [allLabels(end).transition; ...
                find(immTransNoLabel{i},1,'first')];
            allLabels(end).guardempty = [allLabels(end).guardempty; true];
        end
    end

else % here: any(cell2mat(immTransLabel)) = true

    % assume no synchronization label meets necessary criteria
    syncTransitionFound = false;

    % go over immediate transitions with synchronization labels, check if
    % they can be executed -> if yes, find all transition for merge
    for i=1:length(immTransLabel)
        for j=1:length(immTransLabel{i})
    
            % only check immediate transitions
            if immTransLabel{i}(j)
    
                % synchronization label of immediate transition
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
        % no immediate transition with synchronization labels found

        % dummy values for tracker (see if-branch: no immediate transitions)
        tracker(end).transition = zeros(1,length(locID));
        tracker(end).syncLabel = '';

        % quick exit
        restart = false; return
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
% immediate transitions; now get these transitions
for k=1:length(pHA.components)
    if any(allLabels(idxSyncLabel).component == k)
        idxIdReset(k) = false;
    end
end

% initialize transition for full state
transList = cell(1,length(pHA.components));

% loop over all components for full transition
for k=1:length(pHA.components)

    % check for identity reset or not
    if idxIdReset(k)
        % append virtual self-transition with empty guard set and identity
        % reset function
        resetStruct.A = eye(pHA.components(k).location(locID(k)).contDynamics.dim);
        resetStruct.c = zeros(pHA.components(k).location(locID(k)).contDynamics.dim,1);
        transList{k}{1} = transition([],resetStruct,locID(k),syncLabel);

        % target location = same as current location in component k
        list{1}.loc(k) = locID(k);

        % append to list of all labels
        allLabels(idxSyncLabel).component = [allLabels(idxSyncLabel).component; k];
        allLabels(idxSyncLabel).location = [allLabels(idxSyncLabel).location; locID(k)];
        % index this virtual transition by "0" (consistent with livelock check)
        allLabels(idxSyncLabel).transition = [allLabels(idxSyncLabel).transition; 0];
        allLabels(idxSyncLabel).guardempty = [allLabels(idxSyncLabel).guardempty; true];

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
        transList{k}{1} = transition(...
            pHA.components(k).location(locID(k)).transition(idx).guard,...
            pHA.components(k).location(locID(k)).transition(idx).reset,...
            pHA.components(k).location(locID(k)).transition(idx).target,...
            syncLabel);

        % target location
        list{1}.loc(k) = pHA.components(k).location(locID(k)).transition(idx).target;

        % update tracker: idx of transition in transition list of location
        tracker(end).transition = [tracker(end).transition; idx];

    end
    
end


% construct transition for full state (note that the synchronization label
% is removed during the process)
mergedTrans = mergeTransitionSets(pHA,transList,locID,allLabels);

% reset state
list{1}.set = reset(mergedTrans{1},fullR0,fullU);

% other immediate transitions that are now available will be handled in the
% next iteration of the main loop
restart = true;

% notify user that an immediate transition has occurred
if verbose
    disp("  transition: locations [" + ...
        strjoin(string(locID),",") + "] -> locations [" + ...
        strjoin(string(list{1}.loc),",") + "]... " + ...
        "(time: " + string(list{1}.time) + ")");
end

end


% Auxiliary Functions -----------------------------------------------------

function skip = aux_skipTransition(currComp,locID,syncLabel,allLabels)
% check whether current transition is an immediate transition: we have to
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

        % check if the corresponding transition is immediate
        locIdx = compIdx & ...
            allLabels(idxSyncLabel).location == locID(comp);
        if isempty(locIdx) || ~allLabels(idxSyncLabel).guardempty(locIdx)
            skip = true; break;
        end
    end

end

%------------- END OF CODE --------------
