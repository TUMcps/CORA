function han = plotOverTime(simRes,varargin)
% plotOverTime - plots the simulated trajectories over time
%
% Syntax:
%    han = plotOverTime(simRes)
%    han = plotOverTime(simRes,dims)
%    han = plotOverTime(simRes,dims,type)
%
% Inputs:
%    simRes - simResult object
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
% See also: simResult, simulateRandom, plot

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       06-June-2020
% Last update:   22-March-2023 (TL, fixed bug for class arrays)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values for the optional input arguments
dims = setDefaultValues({1},varargin);

% check input arguments
inputArgsCheck({{simRes,'att','simResult','nonempty'};
                {dims,'att','numeric',{'positive','vector','integer'}}});

% check hold status
holdStatus = ishold;
if ~holdStatus
    plot(NaN,NaN,'HandleVisibility','off');
    % reset color index (before readPlotOptions!)
    set(gca(),'ColorOrderIndex',1);
end

% parse input arguments
NVpairs = readPlotOptions(varargin(2:end),'simResult');
[NVpairs,whichtraj] = readNameValuePair(NVpairs,'Traj','ischar','x');

% check which trajectory has to be plotted
whichtraj = aux_checkTraj(simRes,whichtraj);

% save color index
ax = gca();
oldColorIndex = ax.ColorOrderIndex;

% loop over all simulated trajectories
hold on
for i = 1:length(simRes)
    sim_i = simRes(i);
    for j = 1:length(sim_i.(whichtraj))
        han_i = plot(sim_i.t{j},sim_i.(whichtraj){j}(:,dims),NVpairs{:});
    
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

function whichtraj = aux_checkTraj(simRes,whichtraj)

% must be character vector for switch-expression to work properly
if isempty(whichtraj)
    % default value
    whichtraj = 'x';
end

switch whichtraj
    case 'x'
        % no issues (should always be there)
        
    case 'y'
        if isempty(simRes(1).y)
            throw(CORAerror('CORA:emptyProperty'));
        end

    case 'a'
        if isempty(simRes(1).a)
            throw(CORAerror('CORA:emptyProperty'));
        end

    otherwise
        % error
        throw(CORAerror('CORA:specialError','Wrong value for name-value pair'));
end

end

% ------------------------------ END OF CODE ------------------------------
