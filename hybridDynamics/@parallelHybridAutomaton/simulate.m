function [t,x,loc] = simulate(pHA,params)
% simulate - simulates a parallel hybrid automaton
%
% Syntax:
%    [t,x,loc] = simulate(pHA,params)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    params - system parameters
%
% Outputs:
%    t - cell-array storing the time vectors
%    x - cell-array storing the state trajectories
%    loc - double-array storing the visited locations
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: simulateRandom

% Authors:       Victor Charlent, Johann Schoepfer, Niklas Kochdumper
% Written:       24-May-2016  
% Last update:   04-July-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % new options preprocessing
    options = validateOptions(pHA,mfilename,params,struct([]));

    % initialization
    tInter = options.tStart;         % intermediate time at transitions
    locCurr = options.startLoc;      % current location
    xInter = options.x0;             % intermediate state at transitions

    t = {}; x = {}; loc = [];
    cnt = 0;

    % create list of label occurences to check whether all labeled
    % transitions are enabled at the same time
    allLabels = labelOccurrences(pHA);

    % loop over the different locations 
    while tInter < options.tFinal && ~isempty(locCurr) && ...
           ~all(options.finalLoc == locCurr)

        % increment counter
        cnt = cnt + 1;

        % construct new location with local Automaton Product
        currentLocation = locationProduct(pHA,locCurr,allLabels);

        % construct system input for this location
        params.u = aux_mergeInputVector(locCurr,options);

        % simulate within the current location
        params.tStart = tInter;
        params.x0 = xInter;
        params.loc = locCurr;

        [tNew,xNew,locCurr,xInter] = simulate(currentLocation,params);

        % update time
        tInter = tNew(end);

        % concatenate to one full trajectory via cell-arrays (we cannot
        % splice the individual parts into one large sequence since the
        % number of states can vary across locations)
        t{cnt,1} = tNew;
        x{cnt,1} = xNew;
        loc(cnt,:) = params.loc';
    end
    
end


% Auxiliary functions -----------------------------------------------------

function res = aux_mergeInputVector(loc,options)
% construct the input vector for the current location from the inputs for 
% the subcomponents

    res = zeros(size(options.inputCompMap,1),1);
    temp = unique(options.inputCompMap);
    
    for i = 1:length(temp)
        res(options.inputCompMap == temp(i)) = ...
            options.uLoc{temp(i)}{loc(temp(i))}; 
    end
end

% ------------------------------ END OF CODE ------------------------------
