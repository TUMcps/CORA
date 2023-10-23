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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% new options preprocessing
params = validateOptions(HA,mfilename,params,struct([]));

% initialization
tInter = params.tStart;         % intermediate time at transitions
locCurr = params.startLoc;      % current location
xInter = params.x0;             % intermediate state at transitions

% init output arguments
t = {}; x = {}; loc = [];
cnt = 0;

% iteratively simulate for each location separately
while tInter < params.tFinal && ...
      ~isempty(locCurr) && ~isFinalLocation(locCurr,params.finalLoc)

    % increment counter
    cnt = cnt + 1;

    % choose input
    params.u = params.uLoc{locCurr};

    % simulate within the current location
    params.tStart = tInter;
    params.x0 = xInter;
    params.loc = locCurr;
    [tNew,xNew,locCurr,xInter] = simulate(HA.location(locCurr),params);

    % update time
    tInter = tNew(end);

    % concatenate to one full trajectory via cell-arrays (we cannot splice
    % the individual parts into one large sequence since the number of
    % states can vary across locations)
    t{cnt,1} = tNew;
    x{cnt,1} = xNew;
    loc(cnt,1) = params.loc;
end

% ------------------------------ END OF CODE ------------------------------
