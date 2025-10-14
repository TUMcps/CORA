function han = plotTimeStep(traj,varargin)
% plotTimeStep - plots the time step size used in individual
%    simulations over time (all simulations in one graph)
%
% Syntax:
%    han = plotTimeStep(traj)
%    han = plotTimeStep(traj,type)
%
% Inputs:
%    traj - trajectory object
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

% Authors:       Mark Wetzlinger, Laura Luetzow
% Written:       08-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{traj,'att','trajectory','nonempty'}});

% check hold status
holdStatus = ishold;
if ~holdStatus
    plot(NaN,NaN,'HandleVisibility','off');
    % reset color index (before readPlotOptions!)
    set(gca(),'ColorOrderIndex',1);
end

% parse plot options
NVpairs = readPlotOptions(varargin(1:end),'trajectory');

% min / max for axis (if time is const, eps differences are shown...)
mintimestep = Inf; maxtimestep = -Inf;
cumsummin = Inf; cumsummax = -Inf;

% save color index
ax = gca();
oldColorIndex = ax.ColorOrderIndex;

% loop over all simulations
nrSim = length(traj);
hold on; box on;
for r=1:nrSim
    % time axis
    cumsumtVec = [traj(r).t(1,1) repelem(traj(r).t(1,2:end-1),2) traj(r).t(1,end)];
    % time step sizes
    tVec = repelem(diff(traj(r).t(1,:)),2);
    % plot
    plot(cumsumtVec,tVec,NVpairs{:});

    % for axis limits
    if min(tVec) < mintimestep;     mintimestep = min(tVec);     end
    if max(tVec) > maxtimestep;     maxtimestep = max(tVec);     end
    if cumsumtVec(1) < cumsummin;   cumsummin = cumsumtVec(1);   end
    if cumsumtVec(end) > cumsummax; cumsummax = cumsumtVec(end); end
end

% correct color index
updateColorIndex(oldColorIndex);

% labels
xlabel('$t$','interpreter','latex');
ylabel('$\Delta t$','interpreter','latex');
% axes
axis([cumsummin,cumsummax,0.9*mintimestep,1.1*maxtimestep]);

% get handle for graphics object
han = get(groot,'CurrentFigure');

% reset hold status
if ~holdStatus
    hold off
end

if nargout == 0
    clear han;
end

% ------------------------------ END OF CODE ------------------------------
