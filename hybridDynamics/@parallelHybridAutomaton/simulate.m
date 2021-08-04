function [t,x,loc] = simulate(obj,params)
% simulate - simulates a parallel hybrid automaton
%
% Syntax:  
%    [t,x,loc] = simulate(obj, options)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    params - system parameters
%    options - simulation options
%
% Outputs:
%    t - cell-array storing the time vectors
%    x - cell-array storing the state trajectories
%    loc - cell-array storing the visited locations
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: simulateRandom

% Author:        Victor Charlent, Johann Schoepfer, Niklas Kochdumper
% Written:       24-May-2016  
% Last update:   04-July-2018
% Last revision: ---

%------------- BEGIN CODE --------------

    % new options preprocessing
    options = validateOptions(obj,mfilename,params,struct([]));

    % initialization
    tInter = options.tStart;         % intermediate time at transitions
    locCurr = options.startLoc;      % current location
    xInter = options.x0;             % intermediate state at transitions

    loc = {}; t = {}; x = {};

    % loop over the different locations 
    while tInter < options.tFinal && ~isempty(locCurr) && ...
           ~all(options.finalLoc == locCurr)

        % construct new location with local Automaton Product
        currentLocation = locationProduct(obj,locCurr);

        % construct system input for this location
        params.u = mergeInputVector(locCurr,options);

        % simulate within the current location
        params.tStart = tInter;
        params.x0 = xInter;
        params.loc = locCurr;

        [tNew,xNew,locCurr,xInter] = simulate(currentLocation,params);

        % update time
        tInter = tNew(end);

        % store results
        t{end+1,1} = tNew; x{end+1,1} = xNew; loc{end+1,1} = params.loc;
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = mergeInputVector(loc,options)
% construct the input vector for the current location from the inputs for 
% the subcomponents

    res = zeros(size(options.inputCompMap,1),1);
    temp = unique(options.inputCompMap);
    
    for i = 1:length(temp)
       res(options.inputCompMap == temp(i)) = options.uLoc{temp(i)}{loc(temp(i))}; 
    end
end

%------------- END OF CODE --------------