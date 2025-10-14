function han = plotOverTime(traj,varargin)
% plotOverTime - plots the simulated trajectories over time
%
% Syntax:
%    han = plotOverTime(traj)
%    han = plotOverTime(traj,dims)
%    han = plotOverTime(traj,dims,type)
%
% Inputs:
%    traj - trajectory object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%          'Traj', <whichtraj> corresponding to
%                   x ... state trajectory (default)
%                   y ... output trajectory (default if no ti)
%                   a ... algebraic trajectory
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: trajectory, simulateRandom, plot

% Authors:       Niklas Kochdumper, Tobias Ladner, Laura Luetzow
% Written:       08-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values for the optional input arguments
dims = setDefaultValues({1},varargin);

% check input arguments
inputArgsCheck({{traj,'att','trajectory','nonempty'};
                {dims,'att','numeric',{'positive','vector','integer'}}});

% check hold status
holdStatus = ishold;
if ~holdStatus
    plot(NaN,NaN,'HandleVisibility','off');
    % reset color index (before readPlotOptions!)
    set(gca(),'ColorOrderIndex',1);
end

% parse input arguments
NVpairs = readPlotOptions(varargin(2:end),'trajectory');
[NVpairs,whichtraj] = readNameValuePair(NVpairs,'Traj','ischar','x');

% check which trajectory has to be plotted
whichtraj = aux_checkTraj(traj,whichtraj);

% save color index
ax = gca();
oldColorIndex = ax.ColorOrderIndex;

% loop over all simulated trajectories
hold on
for i = 1:length(traj)
    traj_i = traj(i);
    for j = 1:size(traj_i.(whichtraj),3)
        han_i = plot(traj_i.t(1,:),traj_i.(whichtraj)(dims,:,j),NVpairs{:});
    
        if i == 1 && j == 1
            han = han_i;
            % don't display subsequent plots in legend
            NVpairs = [NVpairs, {'HandleVisibility','off'}];
        end
    end
end

% correct color index
updateColorIndex(oldColorIndex);

% reset hold status
if ~holdStatus
    hold off
end

if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function whichtraj = aux_checkTraj(traj,whichtraj)

% must be character vector for switch-expression to work properly
if isempty(whichtraj)
    % default value
    whichtraj = 'x';
end

switch whichtraj
    case 'u'
        if isempty(traj(1).u)
            throw(CORAerror('CORA:emptyProperty'));
        end

    case 'x'
        if isempty(traj(1).x)
            throw(CORAerror('CORA:emptyProperty'));
        end

    case 'y'
        if isempty(traj(1).y)
            throw(CORAerror('CORA:emptyProperty'));
        end

    case 'a'
        if isempty(traj(1).a)
            throw(CORAerror('CORA:emptyProperty'));
        end

    otherwise
        % error
        throw(CORAerror('CORA:specialError','Wrong value for name-value pair'));
end

end

% ------------------------------ END OF CODE ------------------------------
