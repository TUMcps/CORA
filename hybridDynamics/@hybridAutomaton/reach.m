function [R,res] = reach(HA,params,options,varargin)
% reach - computes the reachable set of a hybrid automaton
%
% Syntax:
%    R = reach(HA,params,options)
%    [R,res] = reach(HA,params,options,spec)
%
% Inputs:
%    HA - hybridAutomaton object
%    params - parameter defining the reachability problem
%    options - options for the computation of the reachable set
%    spec - (optional) object of class specification
%
% Outputs:
%    R - reachSet object storing the reachable set
%    res - true/false whether specifications are satisfied
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location/reach

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       07-May-2007 
% Last update:   16-August-2007
%                20-August-2013
%                30-October-2015
%                22-August-2016
%                19-December-2019 (NK, restructured the algorithm)
%                13-October-2021 (MP, location-specific specifications)
%                27-November-2022 (MW, restructure specification syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    res = true;
    
    % options preprocessing
    options = validateOptions(HA,mfilename,params,options);
    % ensure that options are not checked again (in any contDynamics/reach)
    options.validated = true;

    % compute derivatives for each location
    compDerivatives(HA,options);
    
    % check specifications
    spec = [];
    if nargin >= 4
       spec = varargin{1};
    end
    options = aux_check_flatHA_specification(options,HA,spec);

    % initialize reachable set: we use a fixed size to start with and then
    % double the size if the current size is exceeded; this process avoids
    % costly 'and-1' concatenation
    R(10,1) = reachSet();
    % index to append new reachSet objects to full list
    r = 0;

    % initialize queue for reachable set computation (during the
    % computation, multiple branches of reachable set can emerge, requiring
    % to compute the successor reachable sets for all branches one after
    % the other; the queue handles this process)
    list.set = options.R0;
    list.loc = options.startLoc;
    list.time = interval(options.tStart);
    list.parent = 0;

    % display information on command window
    if options.verbose
        disp(newline + "Start analysis...");
    end

    % loop until the queue is empty or a specification is violated
    while ~isempty(list) && res

        % get location, initial set, start time, and parent branch for
        % reachable set computation of first element in the queue
        locID = list(1).loc;
        R0 = list(1).set;
        tStart = list(1).time;
        parent = list(1).parent;

        % get inputs and specification for the current location
        options.U = options.Uloc{locID};
        options.u = options.uloc{locID};
        options.specification = options.specificationLoc{locID};
        % get timeStep for the current location (unless adaptive)
        if ~strcmp(options.linAlg,'adaptive')
            options.timeStep = options.timeStepLoc{locID};
        end

        % check if current location has an instant transition
        instantTransition = arrayfun(@(x) isa(x.guard,'fullspace'),...
            HA.location(locID).transition,'UniformOutput',true);

        if any(instantTransition)

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

            % append to the end of the list
            list(1).set = reset(HA.location(locID).transition(instantTransition),...
                            R0,options.U);
            list(1).loc = HA.location(locID).transition(instantTransition).target;
            list(1).parent = r;

            % print on command window that an instant transition has occurred
            if options.verbose
                disp("  transition: location " + locID + ...
                    " -> location " + list(1).loc + "..." + ...
                    " (time: " + string(tStart) + ")");
            end
            continue

        else
            if options.verbose
                fprintf("Compute reachable set in location " + locID + ...
                    "..." + " (time: " + string(tStart) + " to ");
            end

            % compute the reachable set within a location until either the
            % final time is reached or the reachable set hits a guard set
            % and the computation proceeds in another location
            [Rtemp,Rjump,res] = reach(HA.location(locID),R0,tStart,options);

            if options.verbose
                fprintf(string(Rtemp(1).timePoint.time{end}) + ")\n");
            end
        end

        % remove current element from the queue
        list = list(2:end);

        % add the new branches of reachable sets to the queue
        for i=1:length(Rjump)
            Rjump(i).parent = Rjump(i).parent + r;
        end
        list = [list; Rjump];

        % display transitions on command window
        if options.verbose
            aux_displayTransition(Rjump,locID);
        end

        % store the computed reachable set
        for i=1:size(Rtemp,1)
            % compute output set
            Ytemp = aux_outputSet(Rtemp(i),HA.location(locID),options);
            % init reachSet object and append to full list
            temp = reachSet(Ytemp.timePoint,Ytemp.timeInterval,...
                            parent,locID);
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

    if options.verbose
        disp("...time horizon reached, analysis finished." + newline);
    end
end


% Auxiliary functions -----------------------------------------------------

function options = aux_check_flatHA_specification(options,HA,spec)
% rewrites specifications in the correct format

    numLoc = length(HA.location);

    % initialize specifications with empty cells
    if isempty(spec)

        options.specificationLoc = cell(numLoc,1);

    % adjust specification for each location  
    elseif ~iscell(spec)

        % number of specifications
        nrSpecs = length(spec);

        % checks
        for i=1:length(spec)
            % ensure that time information is not provided (unsupported)
            if ~representsa_(spec(i).time,'emptySet',eps)
                throw(CORAerror('CORA:notSupported',...
                    'Timed specifications are not yet supported for hybrid automata!')); 
            end
            % ensure that no specification is active in a non-existing location
            if any(spec(i).location > numLoc)
                throw(CORAerror('CORA:wrongValue','fourth',...
                    'spec.location must not exceed the number of locations in the hybrid automaton.')); 
            end
        end

        
        options.specificationLoc = cell(numLoc,1);
        % if spec.location = [], specification is assumed to be active in
        % all locations

        for i=1:numLoc
            for j=1:nrSpecs
                if isemptyobject(spec(j).location) || any(spec(j).location == i)
                    options.specificationLoc{i} = ...
                        add(options.specificationLoc{i},spec(j));
                end
            end
        end

    % copy specification for each location
    else

        % check if time information is provided
        for i = 1:length(spec)
            for j = 1:length(spec{i})
                if ~representsa_(spec{i}(j).time,'emptySet',eps)
                    throw(CORAerror('CORA:notSupported',...
                        'Timed specifications are not yet supported for hybrid automata!')); 
                end
            end
        end
        
        % copy specifications
        if all(size(spec) ~= [numLoc,1])
            throw(CORAerror('CORA:notSupported',...
                'Input argument "spec" has the wrong format!'));
        end

        options.specificationLoc = spec;
    end
end

function Ytemp = aux_outputSet(Rtemp,loc,options)
% since we require the reachable set in the entire computation due to guard
% intersections and the preparation of the start set for the next location,
% we only compute the output set at the end of the analysis of each
% location (using the outputSet-functions in contDynamics)

% rewrite options.u of current location to options.uTrans for outputSet()
% TODO: do this in validateOptions?
if isfield(options,'u')
    options.uTrans = options.u;
end

% init
Ytemp.timePoint = [];
Ytemp.timeInterval = [];

% time-point solution
if ~isempty(Rtemp.timePoint)    
    % time-point solution
    nrTimePointSets = length(Rtemp.timePoint.set);
    Ytemp.timePoint.set = cell(nrTimePointSets,1);
    Ytemp.timePoint.time = Rtemp.timePoint.time;
    for i=1:length(Rtemp.timePoint.set)
        Ytemp.timePoint.set{i} = ...
            outputSet(loc.contDynamics,options,Rtemp.timePoint.set{i});
    end
end

% time-interval solution
if ~isempty(Rtemp.timeInterval)
    nrTimeIntervalSets = length(Rtemp.timeInterval.set);
    Ytemp.timeInterval.set = cell(nrTimeIntervalSets,1);
    Ytemp.timeInterval.time = Rtemp.timeInterval.time;
    for i=1:length(Rtemp.timeInterval.set)
        Ytemp.timeInterval.set{i} = ...
            outputSet(loc.contDynamics,options,Rtemp.timeInterval.set{i});
    end
end

% parent
Ytemp.parent = Rtemp.parent;

end

function aux_displayTransition(list,locID)
% only if verbose = true: print outgoing transitions with target location
% identifier and time during which the transition has occurred

if isempty(list)
    return

elseif isscalar(list)
    disp("  transition: location " + strjoin(string(locID),",") + ...
        " -> location " + strjoin(string(list(1).loc),",") + ...
        "... (time: " + string(list(1).time) + ")");

else
    % print 'header'
    fprintf("  transitions: ");
    % indent not in first line
    indent = "";

    % loop over multiple transitions
    for i=1:length(list)
        fprintf(indent + "location " + strjoin(string(locID),",") + ...
            " -> location " + strjoin(string(list(i).loc),",") + ...
            "... (time: " + string(list(i).time) + ")\n");
        % add indent for vertical alignment to all other transitions
        indent = "               ";
    end

end

end

% ------------------------------ END OF CODE ------------------------------
