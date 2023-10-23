function han = plotTimeStep(R,varargin)
% plotTimeStepVector - plots the time step size used in reachSet object
%    over time (all in one graph)
%
% Syntax:
%    han = plotTimeStep(R)
%    han = plotTimeStep(R,type)
%
% Inputs:
%    R - reachSet object
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: @simResult/plotTimeStep

% Authors:       Mark Wetzlinger
% Written:       18-June-2020
% Last update:   15-February-2022 (add case where no R.timePoint given)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse plot options
NVpairs = readPlotOptions(varargin(2:end));

% min / max for axis (if time is const, eps differences are shown...)
mintimestep = Inf; maxtimestep = -Inf; cumsummin = Inf; cumsummax = -Inf;

% save color index
ax = gca();
oldColorIndex = ax.ColorOrderIndex;

hold on; box on;
% loop over all branches in R
for i=1:size(R,1)
    if ~isempty(R(i).timePoint)
        % time axis
        cumsumtVec = [0;repelem(cell2mat(R(i).timePoint.time(1:end-1)),2);...
            R(i).timePoint.time{end}];
        % time step sizes
        tVec = repelem(diff([0;cell2mat(R(i).timePoint.time)]),2);
    else % use timeInterval
        tVec = zeros(length(R(i).timeInterval.time),1);
        for j=1:length(R(i).timeInterval.time)
            tVec(j) = 2*rad(R(i).timeInterval.time{j});
        end
        cumsumtVec = [0;repelem(cumsum(tVec(1:end-1)),2);sum(tVec)];
        tVec = repelem(tVec,2);
    end
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

% title and labels
title('ReachSet: Time Step Size');
xlabel('t');
ylabel('Time Step Size');
% axes
axis([cumsummin,cumsummax,0.9*mintimestep,1.1*maxtimestep]);

% get handle for graphics object
han = get(groot,'CurrentFigure');

if nargout == 0
    clear han;
end

% ------------------------------ END OF CODE ------------------------------
