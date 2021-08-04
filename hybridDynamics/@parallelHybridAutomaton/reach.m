function [R,res] = reach(obj,params,options,varargin)
% reach - computes the reachable set for a parallel hybrid automaton
%
% Syntax:  
%    R = reach(obj,params,options)
%    [R,res] = reach(obj,params,options,spec)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    params - parameter defining the reachability problem
%    options - options for the computation of reachable sets
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
% See also: parallelHybridAutomaton

% Author:       Niklas Kochdumper
% Written:      04-July-2018 
% Last update:  14-June-2020
% Last revision: ---

%------------- BEGIN CODE --------------

    res = 1;
    
    % new options preprocessing
    options = validateOptions(obj,mfilename,params,options);

    options.specification = [];
    if nargin >= 4
       options.specification = varargin{1}; 
    end
    
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
        
        % construct new location with local Automaton Product
        locObj = locationProduct(obj,locID);

        % get input, time step, and specification for the current location
        options.U = mergeInputSet(locID,options);

        % compute the reachable set within a location
        [Rtemp,Rjump,res] = reach(locObj,R0,tStart,options);

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
    
    
% Auxiliary Functions -----------------------------------------------------

function U = mergeInputSet(loc,options)

    numInp = length(options.inputCompMap);
    
    % loop over all used components
    comp = unique(options.inputCompMap);
    
    for i = 1:length(comp)
        
       % find indizes of component in input set
       ind = find(options.inputCompMap == comp(i));
       
       % project to higher dimensional space
       Utemp = projectHighDim(options.Uloc{comp(i)}{loc(comp(i))},numInp,ind);
       
       % add to overall input set
       if i == 1
          U = Utemp; 
       else
          U = U + Utemp; 
       end
    end
end

%------------- END OF CODE --------------