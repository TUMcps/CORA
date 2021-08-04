function [R,res] = reach(obj,params,options,varargin)
% reach - computes the reachable set of a hybrid automaton
%
% Syntax:  
%    R = reach(obj,params,options)
%    [R,res] = reach(obj,params,options,spec)
%
% Inputs:
%    obj - hybrid automaton object
%    params - parameter defining the reachability problem
%    options - options for the computation of the reachable set
%    spec - object of class specification 
%
% Outputs:
%    R - object of class reachSet storing the reachable set
%    res  - 1 if specifications are satisfied, 0 if not
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
% Last revision: ---

%------------- BEGIN CODE --------------

    res = 1;
    
    % options preprocessing
    options = validateOptions(obj,mfilename,params,options);

    % compute derivatives for each location
    compDerivatives(obj,options);
    
    spec = [];
    if nargin >= 4
       spec = varargin{1}; 
    end
    options = check_flatHA_specification(options,obj,spec);
    
    % initialization
    R = [];

    list{1}.set = options.R0;
    list{1}.loc = options.startLoc;
    list{1}.time = interval(options.tStart);
    list{1}.parent = 0;

    % loop until the queue is empty
    while ~isempty(list) && res

        % get first element of the queue
        locID = list{1}.loc;
        R0 = list{1}.set;
        tStart = list{1}.time;
        parent = list{1}.parent;

        list = list(2:end);

        % get input and specification for the current location
        options.U = options.Uloc{locID};
        options.specification = options.specificationLoc{locID};
        % get timeStep for the current location (unless adaptive)
        if ~strcmp(options.linAlg,'adaptive')
            options.timeStep = options.timeStepLoc{locID};
        end

        % compute the reachable set within a location
        [Rtemp,Rjump,res] = reach(obj.location{locID},R0,tStart,options);

        % add the new sets to the queue
        for i = 1:length(Rjump)
            Rjump{i}.parent = Rjump{i}.parent + length(R);
        end
        
        list = [list; Rjump];

        % store the computed reachable set
        for i = 1:size(Rtemp,1)
           temp = reachSet(Rtemp.timePoint,Rtemp.timeInterval, ...
                           Rtemp.parent,locID);
           R = add(R,temp,parent);
        end
    end
end



% Auxiliary Function ------------------------------------------------------

function options = check_flatHA_specification(options,obj,spec)
% rewrites specifications in the correct format

    locations = obj.location;
    numLoc = length(locations);

    % initialize specifications with empty cells
    if isempty(spec)

        specLoc = cell(numLoc,1);

    % same specification for each location  
    elseif ~iscell(spec)   

        specLoc = cell(numLoc,1);

        for i = 1:numLoc
            specLoc{i} = spec;
        end

    % copy specification for each location
    else

        if all(size(spec) ~= [numLoc,1])
            error('Input argument "spec" has the wrong format!');
        else
            specLoc = spec;
        end
    end

    options.specificationLoc = specLoc;

end

%------------- END OF CODE --------------