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

% Authors:       Mark Wetzlinger
% Written:       18-June-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{simRes,'att','simResult','nonempty'}});

% check hold status
holdStatus = ishold;
if ~holdStatus
    plot(NaN,NaN,'HandleVisibility','off');
    % reset color index (before readPlotOptions!)
    set(gca(),'ColorOrderIndex',1);
end

% parse plot options
NVpairs = readPlotOptions(varargin(1:end));

% min / max for axis (if time is const, eps differences are shown...)
mintimestep = Inf; maxtimestep = -Inf;
cumsummin = Inf; cumsummax = -Inf;

% save color index
ax = gca();
oldColorIndex = ax.ColorOrderIndex;

% loop over all simulations
nrSim = length(simRes);
hold on; box on;
for r=1:nrSim
    for i=1:length(simRes(r))
        % time axis
        cumsumtVec = [simRes(r).t{i}(1);repelem(simRes(r).t{i}(2:end-1),2);simRes(r).t{i}(end)];
        % time step sizes
        tVec = repelem(diff(simRes(r).t{i}),2);
        % plot
        plot(cumsumtVec,tVec,NVpairs{:});
    
        % for axis limits
        if min(tVec) < mintimestep;     mintimestep = min(tVec);     end
        if max(tVec) > maxtimestep;     maxtimestep = max(tVec);     end
        if cumsumtVec(1) < cumsummin;   cumsummin = cumsumtVec(1);   end
        if cumsumtVec(end) > cumsummax; cumsummax = cumsumtVec(end); end
    end
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
