function [t,x,loc] = simulate(HA,params)
% simulate - simulates a trajectory of a hybrid automaton
%
% Syntax:
%    [t,x] = simulate(HA,params)
%    [t,x,loc] = simulate(HA,params)
%
% Inputs:
%    HA - hybridAutomaton object
%    params - system parameters
%
% Outputs:
%    t - time vector
%    x - state trajectories
%    loc - visited locations
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location/simulate

% Authors:       Matthias Althoff
% Written:       03-May-2007 
% Last update:   08-May-2020 (MW, update interface)
%                21-May-2023 (MW, output single trajectory)
%                10-May-2024 (NK, allow time-varying inputs)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% new options preprocessing
params = validateOptions(HA,params,struct());

% catch the case with time-varying inputs
if size(params.uLoc{1},2) > 1
    [x,t,loc] = aux_simulateTimeVaryingInput(HA,params);
    return;
end

% initialization
tInter = params.tStart;         % intermediate time at transitions
locCurr = params.startLoc;      % current location
xInter = params.x0;             % intermediate state at transitions

% init output arguments
t = []; x = []; loc = [];

% iteratively simulate for each location separately
while tInter < params.tFinal && ...
      ~isempty(locCurr) && ~priv_isFinalLocation(locCurr,params.finalLoc)

    % choose input
    params.u = params.uLoc{locCurr};

    % compute derivatives of reset functions for current location
    HA = derivatives(HA,locCurr);

    % simulate within the current location
    params.tStart = tInter;
    params.x0 = xInter;
    params.loc = locCurr;
    [tNew,xNew,locCurr,xInter] = simulate(HA.location(locCurr),params);

    % update time
    tInter = tNew(end);

    % concatenate to one full trajectory
    t = [t tNew];
    if size(x,1) < size(xNew,1)
        % old location has fewer states: fill with NaN
        x = [x; NaN(size(xNew,1)-size(x,1), size(x,2))];
    elseif size(x,1) > size(xNew,1)
        % new location has fewer states: fill with NaN
        xNew = [xNew; NaN(size(x,1)-size(xNew,1), size(xNew,2))];
    end
    x = [x xNew];
    loc = [loc repmat(params.loc, 1, length(tNew))];
end
end


% Auxiliary functions -----------------------------------------------------

function [x,t,loc] = aux_simulateTimeVaryingInput(HA,params)
% loop over all input changes to simulte a hybrid automaton with a
% time-varying input signal
    
    % initialization
    steps = size(params.uLoc{1},2);
    dt = params.tFinal/steps;
    uLoc = params.uLoc;
    params = rmfield(params,'uLoc');
    t = []; x = []; loc = [];

    % loop over all input changes
    for i = 1:steps

        % update simulate setting for current step
        params.tFinal = params.tStart + dt;

        for j = 1:length(uLoc)
            params.u{j} = uLoc{j}(:,i);
        end

        % simulate the system
        [tNew,xNew,locNew] = simulate(HA,params);

        % concatenate to one full trajectory
        t = [t tNew];
        if size(x,1) < size(xNew,1)
            % old location has fewer states: fill with NaN
            x = [x; NaN(size(xNew,1)-size(x,1), size(x,2))];
        elseif size(x,1) > size(xNew,1)
            % new location has fewer states: fill with NaN
            xNew = [xNew; NaN(size(x,1)-size(xNew,1), size(xNew,2))];
        end
        x = [x xNew];
        loc = [loc repmat(locNew, 1, length(tNew))];

        % update simulation settings
        params.tStart = params.tFinal;
        params.x0 = x(:,end);
        params.startLoc = loc(:,end);
    end
end

% ------------------------------ END OF CODE ------------------------------
