function h = plotTimeStep(R,varargin)
% plotTimeStepVector - plots the time step size used in reachSet
%    object over time (all in one graph)
%
% Syntax:  
%    h = plotTimeStep(R)
%    h = plotTimeStep(R,'k',...)
%
% Inputs:
%    R - reachSet object
%    varargin - plotting preferences
%
% Outputs:
%    h - handle for the resulting graphic object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: @simResult/plotTimeStep

% Author:       Mark Wetzlinger
% Written:      18-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% input processing
type{1} = 'b';
% If only more than one argument is passed (plotting preferences)
if nargin > 1
    type = varargin;
end

% loop over all simulations
hold on; box on;
% min / max for axis (if time is const, eps differences are shown...)
mintimestep = Inf; maxtimestep = -Inf; cumsummin = Inf; cumsummax = -Inf;

% loop over all branches in R
for i=1:size(R,1)
    % time axis
    cumsumtVec = [0;repelem(cell2mat(R(i).timePoint.time(1:end-1)),2);...
        R(i).timePoint.time{end}];
    % time step sizes
    tVec = repelem(diff([0;cell2mat(R(i).timePoint.time)]),2);
    % plot
    plot(cumsumtVec,tVec,type{:});
    
    % for axis
    if min(tVec) < mintimestep;     mintimestep = min(tVec);     end
    if max(tVec) > maxtimestep;     maxtimestep = max(tVec);     end
    if cumsumtVec(1) < cumsummin;   cumsummin = cumsumtVec(1);   end
    if cumsumtVec(end) > cumsummax; cumsummax = cumsumtVec(end); end
end

% title and labels
title('ReachSet: Time Step Size');
xlabel('t');
ylabel('Time Step Size');
% axes
axis([cumsummin,cumsummax,0.9*mintimestep,1.1*maxtimestep]);


h = get(groot,'CurrentFigure');

%------------- END OF CODE --------------