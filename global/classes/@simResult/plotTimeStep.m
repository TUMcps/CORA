function han = plotTimeStep(simRes,varargin)
% plotTimeStep - plots the time step size used in individual
%    simulations over time (all simulations in one graph)
%
% Syntax:  
%    han = plotTimeStep(simRes)
%    han = plotTimeStep(simRes,type)
%
% Inputs:
%    simRes - simResult object
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author:       Mark Wetzlinger
% Written:      18-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% parse plot options
NVpairs = readPlotOptions(varargin(2:end));

% min / max for axis (if time is const, eps differences are shown...)
mintimestep = Inf; maxtimestep = -Inf; cumsummin = Inf; cumsummax = -Inf;

% loop over all simulations
nrSim = length(simRes.t);
hold on; box on;
for i=1:nrSim
    % time axis
    cumsumtVec = [simRes.t{i}(1);repelem(simRes.t{i}(2:end-1),2);simRes.t{i}(end)];
    % time step sizes
    tVec = repelem(diff(simRes.t{i}),2);
    % plot
    plot(cumsumtVec,tVec,NVpairs{:});

    % for axis limits
    if min(tVec) < mintimestep;     mintimestep = min(tVec);     end
    if max(tVec) > maxtimestep;     maxtimestep = max(tVec);     end
    if cumsumtVec(1) < cumsummin;   cumsummin = cumsumtVec(1);   end
    if cumsumtVec(end) > cumsummax; cumsummax = cumsumtVec(end); end
end

% title and labels
title('Simulation: Time Step Size');
xlabel('t');
ylabel('Time Step Size');
% axes
axis([cumsummin,cumsummax,0.9*mintimestep,1.1*maxtimestep]);

% get handle for graphics object
han = get(groot,'CurrentFigure');

if nargout == 0
    clear han;
end

%------------- END OF CODE --------------