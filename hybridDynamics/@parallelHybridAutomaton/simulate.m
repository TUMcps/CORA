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

    % options preprocessing
    params = validateOptions(pHA,params,struct());

    % initialization
    tInter = params.tStart;         % intermediate time at transitions
    locCurr = params.startLoc;      % current location
    xInter = params.x0;             % intermediate state at transitions

    t = {}; x = {}; loc = [];
    cnt = 0;

    % create list of label occurences to check whether all labeled
    % transitions are enabled at the same time
    allLabels = labelOccurrences(pHA);

    % loop over the different locations 
    while tInter < params.tFinal && ~isempty(locCurr) && ...
           ~all(params.finalLoc == locCurr)

        % increment counter
        cnt = cnt + 1;

        % construct new location with local Automaton Product
        currentLocation = locationProduct(pHA,locCurr,allLabels);

        % construct system input for this location
        params.u = aux_mergeInputVector(locCurr,params.uLoc,params.inputCompMap);

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

function res = aux_mergeInputVector(loc,uLoc,inputCompMap)
% construct the input vector for the current location from the inputs for 
% the subcomponents

    res = zeros(size(inputCompMap,1),1);
    temp = unique(inputCompMap);
    
    for i = 1:length(temp)
        res(inputCompMap == temp(i)) = uLoc{temp(i)}{loc(temp(i))}; 
    end
end

% ------------------------------ END OF CODE ------------------------------
