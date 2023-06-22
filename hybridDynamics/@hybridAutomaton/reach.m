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
% See also: location/reach

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      07-May-2007 
% Last update:  16-August-2007
%               20-August-2013
%               30-October-2015
%               22-August-2016
%               19-December-2019 (NK, restructured the algorithm)
%               13-October-2021 (MP, implemented location-specific
%                                specifications)
%               27-November-2022 (MW, restructure specification syntax)
% Last revision:---

%------------- BEGIN CODE --------------

    res = true;
    
    % options preprocessing
    options = validateOptions(HA,mfilename,params,options);

    % compute derivatives for each location
    compDerivatives(HA,options);
    
    % check specifications
    spec = [];
    if nargin >= 4
       spec = varargin{1};
    end
    options = aux_check_flatHA_specification(options,HA,spec); %,spec_locations);

    % initialize reachable set
    R = [];

    % initialize queue for reachable set computation (during the
    % computation, multiple branches of reachable set can emerge, requiring
    % to compute the successor reachable sets for all branches one after
    % the other; the queue handles this process)
    list{1}.set = options.R0;
    list{1}.loc = options.startLoc;
    list{1}.time = interval(options.tStart);
    list{1}.parent = 0;

    % display information on command window
    if options.verbose
        disp(newline + "Start analysis...");
    end

    % loop until the queue is empty or a specification is violated
    while ~isempty(list) && res

        % get location, initial set, start time, and parent branch for
        % reachable set computation of first element in the queue
        locID = list{1}.loc;
        R0 = list{1}.set;
        tStart = list{1}.time;
        parent = list{1}.parent;

        % get input and specification for the current location
        options.U = options.Uloc{locID};
        options.specification = options.specificationLoc{locID};
        % get timeStep for the current location (unless adaptive)
        if ~strcmp(options.linAlg,'adaptive')
            options.timeStep = options.timeStepLoc{locID};
        end

        % check if current location has an immediate transition (this is
        % the case if any guard set is defined as guard = []; also, the 
        % transition must not be empty)
        immediateTransition = arrayfun(@(x) isnumeric(x.guard) && isempty(x.guard) && ~isempty(x.target),...
            HA.location(locID).transition,'UniformOutput',true);

        if any(immediateTransition)
            % execute reset function of immediate transition: overwrite
            % start set and location (time and parent do not change)
            list{1}.set = reset(HA.location(locID).transition(immediateTransition),...
                            R0,options.U);
            list{1}.loc = HA.location(locID).transition(immediateTransition).target;

            % notify user that an immediate transition has occurred
            if options.verbose
                disp("  transition: location " + locID + ...
                    " -> location " + list{1}.loc + "..." + ...
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
            Rjump{i}.parent = Rjump{i}.parent + length(R);
        end
        list = [list; Rjump];

        % display transitions on command window
        if options.verbose && isscalar(list)
            % multiple successor locations currently not supported...
            disp("  transition: location " + locID + " -> location " + ...
                string(list{1}.loc) + "... (time: " + string(list{1}.time) + ")");
        end

        % store the computed reachable set
        for i=1:size(Rtemp,1)
            % compute output set
            Ytemp = aux_outputSet(Rtemp(i),HA.location(locID),options);
            % store in reachSet object
            temp = reachSet(Ytemp.timePoint,Ytemp.timeInterval,...
                            Ytemp.parent,locID);
            R = add(R,temp,parent);
        end
    end

    if options.verbose
        disp("...time horizon reached, analysis finished." + newline);
    end
end



% Auxiliary Function ------------------------------------------------------

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
            if ~isempty(spec(i).time)
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
                if ~isempty(spec{i}(j).time)
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

    % adapt options (as in location/reach)
    [params,options_] = adaptOptions(loc,options);
    % dummy value for R0...
    params.R0 = zonotope(zeros(loc.contDynamics.dim,1));
    options_ = validateOptions(loc.contDynamics,'reach',params,options_);
    
    % time-point solution
    nrTimePointSets = length(Rtemp.timePoint.set);
    Ytemp.timePoint.set = cell(nrTimePointSets,1);
    Ytemp.timePoint.time = Rtemp.timePoint.time;
    for i=1:length(Rtemp.timePoint.set)
        Ytemp.timePoint.set{i} = ...
            outputSet(loc.contDynamics,options_,Rtemp.timePoint.set{i});
    end
    
    % time-interval solution
    nrTimeIntervalSets = length(Rtemp.timeInterval.set);
    Ytemp.timeInterval.set = cell(nrTimeIntervalSets,1);
    Ytemp.timeInterval.time = Rtemp.timeInterval.time;
    for i=1:length(Rtemp.timeInterval.set)
        Ytemp.timeInterval.set{i} = ...
            outputSet(loc.contDynamics,options_,Rtemp.timeInterval.set{i});
    end
    
    % parent
    Ytemp.parent = Rtemp.parent;

end

%------------- END OF CODE --------------