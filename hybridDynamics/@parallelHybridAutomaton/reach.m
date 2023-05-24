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

% Author:       Niklas Kochdumper
% Written:      04-July-2018 
% Last update:  14-June-2020
%               13-October-2021 (MP, implemented location-specific
%                                specifications)
%               03-March-2022 (MP, implemented synchronization labels)
%               27-November-2022 (MW, restructure specification syntax)
%               20-January-2023 (MW, save already computed location products)
% Last revision:---

%------------- BEGIN CODE --------------

    res = true;
    
    % new options preprocessing
    options = validateOptions(pHA,mfilename,params,options);

    % preprocessing of specifications
    options.specification = [];
    spec = [];
    specLocSpecific = false;
    
    if nargin >= 4
        spec = varargin{1};
        % init that all specs are active everywhere
        options.specification = spec;
        % check if location-specific specs given 
        for i=1:length(spec)
            if isempty(spec(i).location)
                specLocSpecific = true; break
            end
        end
    end
    
    % initialize reachable set
    R = [];

    % initialize queue
    list{1}.set = options.R0;
    list{1}.loc = options.startLoc;
    list{1}.time = interval(options.tStart);
    list{1}.parent = 0;

    % display information on command window
    if options.verbose
        disp("Start analysis...");
        disp("  pre-process synchronization labels...");
    end

    % initialize tracker (for livelock detection)
    tracker = struct(...
        'switchingTime',interval(),...
        'locID',double.empty(0,length(options.startLoc)),...
        'transition',double.empty(0,length(options.startLoc)),...
        'syncLabel','');

    % create list of label occurences to check whether all labeled
    % transitions are enabled at the same time
    allLabels = labelOccurrences(pHA);

    % number of iterations in the main loop
    k = 0;
    
    % loop until the queue is empty or a specification is violated
    while ~isempty(list) && res

        % increment number of iterations
        k = k + 1;

        % update tracking
        tracker(k,1).switchingTime = list{1}.time;
        tracker(k,1).locID = list{1}.loc;

        % get locations for each hybrid automaton, initial set for 
        % automaton product, start time, and parent branch for reachable
        % set computation of first element in the queue
        locID = list{1}.loc;
        R0 = list{1}.set;
        tStart = list{1}.time;
        parent = list{1}.parent;

        % filter the specifications for the ones relevant to the local
        % automaton product (if no specification is location-specific, the
        % specifications are always the same ones as initialized)
        if specLocSpecific
            options.specification = aux_filterSpecifications(spec,locID);
        end
        
        % compute input set for the constructed location
        options.U = aux_mergeInputSet(locID,options.Uloc,options.inputCompMap);
        
        % check for immediate transitions
        [list,tracker,restart] = immediateTransition(pHA,list,locID,...
            allLabels,R0,options.U,tracker,options.verbose);
        % check for livelock
        if checkLivelock(tracker); break; end
        % restart if immediate transition has occurred
        if restart; continue; end

        % location for evaluation via local Automaton Product        
        % check if location product has been computed before
        locProdIdx = cellfun(@(x) all(x == locID),{pHA.locProd.locID});
        if k > 0 && any(locProdIdx)
            % read location product from saved locations
            locObj = pHA.locProd(locProdIdx).location;
        else
            % compute new location product
            if options.verbose
                disp("  compute location product of locations [" + ...
                   strjoin(string(locID),',') + "]...");
            end
            locObj = locationProduct(pHA,locID,allLabels);
            pHA.locProd = [pHA.locProd;...
                struct('location',locObj,'locID',locID)];
        end

        % compute the reachable set within the constructed location
        if options.verbose
            disp("  compute reachable set in locations [" + ...
               strjoin(string(locID),',') + "]...");
        end
        [Rtemp,Rjump,res] = reach(locObj,R0,tStart,options);

        % remove current element from the queue
        list = list(2:end);

        % add the new sets to the queue
        for i = 1:length(Rjump)
            Rjump{i}.parent = Rjump{i}.parent + length(R);
        end
        list = [list; Rjump];

        % display transitions on command window
        if options.verbose && isscalar(list)
            % multiple successor locations currently not supported...
            disp("  transition: locations [" + strjoin(string(locID),",") + ...
                "] -> locations [" + strjoin(string(list{1}.loc),",") + ...
                "]... (time: " + string(list{1}.time) + ")");
        end

        % store the computed reachable set
        for i = 1:size(Rtemp,1)
            temp = reachSet(Rtemp.timePoint,Rtemp.timeInterval,...
                                Rtemp.parent,locID);
            R = add(R,temp,parent);
        end
    end
end
    
    
% Auxiliary Functions -----------------------------------------------------

function U = aux_mergeInputSet(loc,Uloc,inputCompMap)
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
        Utemp = projectHighDim(Uloc{comp(i)}{loc(comp(i))},numInp,ind);
        
        % add to overall input set
        if i == 1
            U = Utemp; 
        else
            U = U + Utemp; 
        end
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
        if any(cellfun(@(x) isempty(x),spec(i).location,'UniformOutput',true))
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

%------------- END OF CODE --------------
