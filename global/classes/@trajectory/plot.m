function han = plot(traj,varargin)
% plot - plots a projection of the simulated trajectories
%
% Syntax:
%    han = plot(traj)
%    han = plot(traj,dims)
%    han = plot(traj,dims,type)
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
% See also: trajectory, simulateRandom

% Authors:       Niklas Kochdumper, Matthias Althoff, Laura Luetzow
% Written:       08-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default input arguments
dims = setDefaultValues({[1,2]},varargin);

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

% parse plotting options
NVpairs = readPlotOptions(varargin(2:end),'trajectory');

NVpairs = [NVpairs,{'NVPAIRS_VALIDATED',true}]; 
[NVpairs,whichtraj] = readNameValuePair(NVpairs,'Traj','ischar','x');

% plot jumps to other locations in dotted line
NVpairs_jump = [NVpairs(1:end-2), {'LineStyle', ':','NVPAIRS_VALIDATED',true}]; 

% check dimension
if length(dims) < 2
    throw(CORAerror('CORA:plotProperties',2));
elseif length(dims) > 3
    throw(CORAerror('CORA:plotProperties',3));
end

% check which trajectory has to be plotted
whichtraj = aux_checkTraj(traj,whichtraj);

% save color index
ax = gca();
oldColorIndex = ax.ColorOrderIndex;

% loop over all simulated trajectories
hold on
for r=1:length(traj)
    % extract indizes for jumps to different locations
    jumps = [0 find(sum(abs(diff(traj(r).loc,1,2)),1) ~= 0) traj(r).n_k];
    for i_jumpEnd = 2:length(jumps)
        % indizes for trajectory part within the same location
        idz = jumps(i_jumpEnd-1)+1:jumps(i_jumpEnd);

        % indizes for jump between locations
        if i_jumpEnd > 2
            idz_jump = jumps(i_jumpEnd-1):jumps(i_jumpEnd-1)+1;
        else
            idz_jump = [];
        end

        % plot trajectory parts
        for j = 1:size(traj(r).(whichtraj),3)
            han_i = plotPolygon(traj(r).(whichtraj)(dims,idz,j),NVpairs{:});
            if ~isempty(idz_jump)
                % plot jump
                han_i = plotPolygon(traj(r).(whichtraj)(dims,idz_jump,j),NVpairs_jump{:});
            end
        end
    end

    if r == 1
        han = han_i;
        % don't display subsequent plots in legend
        NVpairs = [NVpairs, {'HandleVisibility','off'}];
        NVpairs_jump = [NVpairs_jump, {'HandleVisibility','off'}];
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
