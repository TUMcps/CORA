function han = plot(simRes,varargin)
% plot - plots a projection of the simulated trajectories
%
% Syntax:
%    han = plot(simRes)
%    han = plot(simRes,dims)
%    han = plot(simRes,dims,type)
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
% See also: simResult, simulateRandom

% Authors:       Niklas Kochdumper, Matthias Althoff
% Written:       06-June-2020
% Last update:   28-July-2020 (MA, 3D plots added)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default input arguments
dims = setDefaultValues({[1,2]},varargin);

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

% parse plotting options
NVpairs = readPlotOptions(varargin(2:end),'simResult');
[NVpairs,whichtraj] = readNameValuePair(NVpairs,'Traj','ischar','x');

% check dimension
if length(dims) < 2
    throw(CORAerror('CORA:plotProperties',2));
elseif length(dims) > 3
    throw(CORAerror('CORA:plotProperties',3));
end

% check which trajectory has to be plotted
whichtraj = aux_checkTraj(simRes,whichtraj);

% save color index
ax = gca();
oldColorIndex = ax.ColorOrderIndex;

% loop over all simulated trajectories
hold on
for r=1:length(simRes)
    for i=1:length(simRes(r).(whichtraj))
        han_i = plotPolygon(simRes(r).(whichtraj){i}(:,dims)',NVpairs{:});
        
        if i == 1
            han = han_i;
            % don't display subsequent plots in legend
            NVpairs = [NVpairs, {'HandleVisibility','off'}];
        end
    end
end

% correct color index
updateColorIndex(oldColorIndex);

% show 3 dimensions ('hold on' causes projection to first two dimensions
%   to be shown if figure was not 3D before)
if length(dims) == 3
    view(3);
end

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
