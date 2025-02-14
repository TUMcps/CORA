function [R,res] = reach(pHA,params,options,varargin)
% reach - computes the reachable set for a parallel hybrid automaton
%
% Syntax:
%    R = reach(pHA,params,options)
%    [R,res] = reach(pHA,params,options,spec)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    params - parameter defining the reachability problem
%    options - options for the computation of reachable sets
%    spec - object of class specification 
%
% Outputs:
%    R - reachSet object storing the reachable set
%    res - true/false whether specifications are satisfied
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: parallelHybridAutomaton

% Authors:       Niklas Kochdumper
% Written:       04-July-2018 
% Last update:   14-June-2020
%                13-October-2021 (MP, location-specific specifications)
%                03-March-2022 (MP, implemented synchronization labels)
%                27-November-2022 (MW, restructure specification syntax)
%                20-January-2023 (MW, save already computed location products)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    spec = setDefaultValues({[]},varargin);
    res = true;
    
    % options preprocessing
    [params,options] = validateOptions(pHA,params,options);

    % preprocessing of specifications
    if ~isempty(spec)
        % init that all specs are active everywhere
        options.specification = spec;
        specLocSpecific = aux_checkSpecLocSpecific(spec);
    else
        options.specification = [];
        specLocSpecific = false;
    end
    
    % initialize reachable set: we use a fixed size to start with and then
    % double the size if the current size is exceeded; this process avoids
    % costly 'and-1' concatenation
    R(10,1) = reachSet();
    % index to append new reachSet objects to full list
    r = 0;

    % initialize queue
    list.set = params.R0;
    list.loc = params.startLoc;
    list.time = interval(params.tStart);
    list.parent = 0;

    % display information on command window
    aux_verbose_displayStart(options.verbose);

    % initialize tracker (for livelock detection)
    tracker = struct(...
        'switchingTime',interval.empty(1),...
        'locID',double.empty(0,length(params.startLoc)),...
        'transition',double.empty(0,length(params.startLoc)),...
        'syncLabel','');

    % create list of label occurences to check whether all labeled
    % transitions are enabled at the same time
    allLabels = priv_labelOccurrences(pHA);

    % number of iterations in the main loop
    k = 0;
    
    % loop until the queue is empty or a specification is violated
    while ~isempty(list) && res

        % increment number of iterations
        k = k + 1;

        % update tracking
        tracker(k,1).switchingTime = list(1).time;
        tracker(k,1).locID = list(1).loc;

        % get locations for each hybrid automaton, initial set for 
        % automaton product, start time, and parent branch for reachable
        % set computation of first element in the queue
        locID = list(1).loc;
        R0 = list(1).set;
        tStart = list(1).time;
        parent = list(1).parent;

        % filter the specifications for the ones relevant to the local
        % automaton product (if no specification is location-specific, the
        % specifications are always the same ones as initialized)
        if specLocSpecific
            options.specification = aux_filterSpecifications(spec,locID);
        end
        
        % compute input set for the constructed location
        [params.U,params.u] = aux_mergeInputSet(...
            locID,params.Uloc,params.uloc,params.inputCompMap);
        % compute disturbance set and noise set for the constructed location
        [params.W,params.V] = aux_mergeDistNoiseSet(locID,params.Wloc,params.Vloc);
        
        % check for instant transitions
        [mergedTrans,tracker] = priv_instantTransition(pHA,locID,allLabels,tracker);
        
        if false %% priv_checkLivelock(tracker)
            % check for livelock (currently deactivated)
            break
        elseif ~isempty(mergedTrans)
            % restart if instant transition has occurred

            % save merged transition
            if isempty(pHA.mergedTrans) || ...
                ~any( all([pHA.mergedTrans.locID] == locID,1) ...
                & all([pHA.mergedTrans.transID] == tracker(end).transition,1) )
                pHA.mergedTrans = [pHA.mergedTrans; ...
                    struct('transition',mergedTrans,'locID',locID,...
                    'transID',tracker(end).transition)];
            end

            % save reachable set to array
            temp = struct('set',{{R0}},'time',{{tStart}});
            if r == length(R)
                R(2*r,1) = reachSet();
            end
            R(r+1,1) = reachSet(temp,[],list(1).parent,locID);
            % increment counter
            r = r + 1;

            % here, we overwrite the first entry in list and continue the
            % reachability analysis with this set -- in contrast to below,
            % where the new sets are appended at the end of the list
            
            % reset state
            list(1).set = evaluate(mergedTrans.reset,R0,params.U);
            % reset location ID
            list(1).loc = mergedTrans.target;
            % reset parent
            list(1).parent = r;
            
            % print on command window that an instant transition has occurred
            aux_verbose_displayInstantTransition(options.verbose,list,locID);
            continue
        end

        % location for evaluation via local Automaton Product
        % check if location product has been computed before
        locProdIdx = cellfun(@(x) all(x == locID),{pHA.locProd.locID});
        if k > 0 && any(locProdIdx)
            % read location product from saved locations
            locObj = pHA.locProd(locProdIdx).location;
        else
            % compute new location product
            aux_verbose_displayLocationProduct(options.verbose,locID);
            locObj = locationProduct(pHA,locID,allLabels);
            pHA.locProd = [pHA.locProd;...
                struct('location',locObj,'locID',locID)];
        end

        % compute the reachable set within the constructed location
        aux_verbose_displayProgress(options.verbose,locID);
        params.R0 = R0; params.tStart = tStart;
        [Rtemp,Rjump,res] = reach(locObj,params,options);

        % remove current element from the queue
        list = list(2:end);

        % add the new sets to the queue
        for i = 1:length(Rjump)
            Rjump(i).parent = Rjump(i).parent + r;
        end
        list = [list; Rjump];

        % display transitions on command window
        aux_verbose_displayTransition(options.verbose,Rjump,locID);

        % store the computed reachable set
        for i = 1:size(Rtemp,1)
            % construct new reachable set
            temp = reachSet(Rtemp.timePoint,Rtemp.timeInterval,...
                                parent,locID);
            % append new reachable set to full list
            if r == length(R)
                R(2*r,1) = reachSet();
            end
            R(r+1,1) = temp;
            % increment counter
            r = r + 1;
        end
    end

    % truncate reachable set (empty entries at the end due to
    % pre-allocation of memory)
    R = R(1:r,1);

    aux_verbose_displayEnd(options.verbose);
end


% Auxiliary functions -----------------------------------------------------

function specLocSpecific = aux_checkSpecLocSpecific(spec)

% check if location-specific specs given 
for i=1:length(spec)
    if isempty(spec(i).location)
        specLocSpecific = true; break
    end
end

end

function [U,u] = aux_mergeInputSet(loc,Uloc,uloc,inputCompMap)
% compute the joint input set for the location generated via the automaton
% product of all individual hybrid automata

    % loop over all used components
    numInp = length(inputCompMap);
    comp = unique(inputCompMap);
    
    for i=1:length(comp)
        
        % find indices of component in input set
        ind = find(inputCompMap == comp(i));
        
        % project input set of location from individual hybrid automaton
        % to higher-dimensional space of automaton product
        Utemp = projectHighDim_(Uloc{comp(i)}{loc(comp(i))},numInp,ind);
        % same for trajectory
        utemp = zeros(numInp,1);
        utemp(ind) = uloc{comp(i)}{loc(comp(i))};
        
        % add to overall input set
        if i == 1
            U = Utemp; 
            u = utemp;
        else
            U = U + Utemp; 
            u = u + utemp;
        end
    end

end

function [W,V] = aux_mergeDistNoiseSet(locID,Wloc,Vloc)
% current assumption: disturbances/noises for all locations are independent
% from one another

W = Wloc{1}{locID(1)};
V = Vloc{1}{locID(1)};

for i=2:length(locID)
    W = cartProd(W,Wloc{i}{locID(i)});
    V = cartProd(V,Vloc{i}{locID(i)});
end

end

function activeSpecs = aux_filterSpecifications(spec,locID)
% filter the specifications for the ones relevant to the local automaton
% product: if a specific subcomponent is within a location where a
% specification is relevant, the specification is checked in the next
% reachability step

    % number of components
    nrComp = length(locID);

    % number of specifications
    nrSpecs = length(spec);
    
    % start with empty list, collect indices of relevant specifications
    activeSpecs_idx = false(nrSpecs,1);
    
    % loop over components and specifications to see if they are active in
    % the current set of locations (given location per subcomponent)
    for i=1:nrSpecs
        % active in all locations of any subcomponent?
        if any(cellfun(@(x) isemptyobject(x),spec(i).location,'UniformOutput',true))
            activeSpecs_idx(i) = true; continue
        end
        % loop over subcomponents
        for j=1:nrComp
            if any(spec(i).location{j} == locID(j))
                activeSpecs_idx(i) = true; break
            end
        end
    end
    
    % return specification object containing only active specifications
    activeSpecs = spec(activeSpecs_idx);

end

% logging functions below... (only print if options.verbose = true)

function aux_verbose_displayStart(verbose)
% display start message
if ~verbose; return; end

disp("Start analysis...");
disp("  pre-process synchronization labels...");

end

function aux_verbose_displayEnd(verbose)
% display end message
if ~verbose; return; end

disp("...time horizon reached, analysis finished." + newline);

end

function aux_verbose_displayInstantTransition(verbose,list,locID)
% display if instant transition happened
if ~verbose; return; end

if options.verbose
    disp("  transition: locations [" + ...
        strjoin(string(locID),",") + "] -> locations [" + ...
        strjoin(string(list(1).loc),",") + "]... " + ...
        "(time: " + string(list(1).time) + ")");
end

end

function aux_verbose_displayLocationProduct(verbose,locID)
% display which location production is computed next
if ~verbose; return; end

disp("  compute location product of locations [" + ...
   strjoin(string(locID),',') + "]...");

end

function aux_verbose_displayProgress(verbose,locID)
% display which location's reachable set is computed next
if ~verbose; return; end

disp("  compute reachable set in locations [" + ...
   strjoin(string(locID),',') + "]...");

end

function aux_verbose_displayTransition(verbose,list,locID)
% only if verbose = true: print outgoing transitions with target location
% identifier and time during which the transition has occurred
if ~verbose; return; end

if isempty(list)
    return

elseif isscalar(list)
    disp("  transition: locations [" + strjoin(string(locID),",") + ...
        "] -> locations [" + strjoin(string(list(1).loc),",") + ...
        "]... (time: " + string(list(1).time) + ")");

else
    % print 'header'
    fprintf("  transitions: ");
    % indent not in first line
    indent = "";

    % loop over multiple transitions
    for i=1:length(list)
        fprintf(indent + "locations [" + strjoin(string(locID),",") + ...
            "] -> locations [" + strjoin(string(list(1).loc),",") + ...
            "]... (time: " + string(list(1).time) + ")\n");
        % add indent for vertical alignment to all other transitions
        indent = "               ";
    end

end

end

% ------------------------------ END OF CODE ------------------------------
