function [t,x,loc] = simulate(obj,params)
% simulate - simulates a hybrid automaton
%
% Syntax:  
%    obj = simulate(obj,params)
%
% Inputs:
%    obj - hybrid automaton object
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
% See also: none

% Author:        Matthias Althoff
% Written:       03-May-2007 
% Last update:   08-May-2020 (MW, update interface)
% Last revision: ---

%------------- BEGIN CODE --------------

% new options preprocessing
params = validateOptions(obj,mfilename,params,struct([]));

% initialization
tInter = params.tStart;         % intermediate time at transitions
locCurr = params.startLoc;      % current location
xInter = params.x0;             % intermediate state at transitions

loc = {}; t = {}; x = {};

% iteratively simulate for each location seperately
while (tInter < params.tFinal) && ...
      (~isempty(locCurr)) && ~isFinalLocation(locCurr,params.finalLoc)

    % choose input
    params.u = params.uLoc{locCurr};

    % simulate within the current location
    params.tStart = tInter;
    params.x0 = xInter;
    params.loc = locCurr;
    [tNew,xNew,locCurr,xInter] = simulate(obj.location{locCurr},params);

    % update time
    tInter = tNew(end);

    % store results
    t{end+1,1} = tNew; x{end+1,1} = xNew; loc{end+1,1} = params.loc;
end

%------------- END OF CODE --------------